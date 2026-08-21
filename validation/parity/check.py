"""Check BetterBuilder against captured oracle values, stage by stage.

    validation/parity/check.py [stage1] [stage2] [stage3]      (default: all)

Exit status 0 if every stage passes. Fixtures come from the capture scripts in this directory; see
README.md for how they are produced and why each metric is shaped the way it is.

Requires the release binaries (`bb-spec`, `bb-embed`) and numpy.
"""
import json
import math
import os
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
FIX = HERE / "fixtures"

# The corpus is an input, not an assumption: a gate can only report on the chemistry it was run
# against. Each is a glob so coverage extends by adding a provenanced corpus and its fixture,
# not by editing this file.
CORPUS = sorted((REPO / "validation").glob(os.environ.get("BB_CORPUS", "seeds_*.smi")))
PERCEPTION = FIX / os.environ.get("BB_PERCEPTION", "perception_seeds.json")
_ENS = os.environ.get("BB_ENSEMBLES", "ensembles_seeds")
# The Divya reference ensembles are large + regenerated on demand (validation/fixtures/divya/), so
# allow an absolute path here; otherwise resolve relative to the committed fixtures dir.
ENSEMBLES = pathlib.Path(_ENS) if os.path.isabs(_ENS) else FIX / _ENS
STAGE2_GLOB = os.environ.get("BB_STAGE2", "stage2_*.json")
EXO_GLOB = os.environ.get("BB_EXO", "exo_*.json")
BB_SPEC = REPO / "rust/target/release/bb-spec"
BB_SPEC_NATIVE = REPO / "rust/target/release/bb-spec-native"
BB_EMBED = REPO / "rust/target/release/bb-embed"
BB_SOLVATE = REPO / "rust/target/release/bb-solvate"
BB_PERCEIVE = REPO / "rust/target/release/bb-perceive"
BB_BOUNDS_DUMP = REPO / "rust/target/release/bb-bounds-dump"

BOUNDS_TOL = 2e-3   # f32 storage of a f64 bounds matrix
V_TOL = 1e-3        # torsion force constants


def mol_id(name):
    """Corpus key for a molecule name. Her records carry a `.<prot_id>` suffix the corpus lacks."""
    base = name.split("_")[0]
    head, sep, tail = base.rpartition(".")
    return head if sep and tail.isdigit() else base


def smiles_by_id():
    out = {}
    for path in CORPUS:
        for line in open(path):
            f = line.split()
            if len(f) >= 2:
                out[mol_id(f[1])] = f[0]
    return out


def spec_for(smiles):
    r = subprocess.run([str(BB_SPEC), smiles], stdout=subprocess.PIPE, universal_newlines=True)
    return json.loads(r.stdout)


def native_spec_for(smiles):
    """The shipped RDKit-free spec (bb-spec-native). Identical to the bridge spec on every gated
    field, but it also carries `raw_bounds` — which the core-pin recipe needs to tighten the
    sidechain bounds the way RDKit's coordMap path does. stage3 embeds this, not the bridge spec, so
    it exercises the pipeline we actually ship rather than native-embed on an oracle spec missing a
    native-only field."""
    r = subprocess.run([str(BB_SPEC_NATIVE), smiles], stdout=subprocess.PIPE, universal_newlines=True)
    return json.loads(r.stdout)


def embed(spec, seed):
    r = subprocess.run([str(BB_EMBED), "-", "recipe", str(seed)], input=json.dumps(spec),
                       stdout=subprocess.PIPE, universal_newlines=True)
    return json.loads(r.stdout)


def iter_ensembles(ens_dir, smi=None):
    """Yield (mid, smiles, spec, her, Hf, path) for each `*.confs.json` in `ens_dir` whose molecule id
    is in the corpus. `Hf` is her conformers reshaped `(K, n_atoms, 3)`. This is the single source for
    the load-ensemble → resolve-native-spec → reshape prologue that recurred verbatim across the
    criterion-2 / TFD / strain scripts. Iteration is by matched molecule (unmatched ids skipped);
    callers apply their own downstream skips (atom-count match vs the native spec, minimum conformer
    count) and their own cap on the number of molecules."""
    if smi is None:
        smi = smiles_by_id()
    for path in sorted(pathlib.Path(ens_dir).glob("*.confs.json")):
        her = json.load(open(path))
        mid = mol_id(her["name"])
        if mid not in smi:
            continue
        spec = native_spec_for(smi[mid])
        Hf = np.array(her["conformers"]).reshape(-1, her["n_atoms"], 3)
        yield mid, smi[mid], spec, her, Hf, path


def heavy(z):
    return [i for i, x in enumerate(z) if x > 1]


def automorphisms(smiles, n_heavy):
    """Heavy-atom graph automorphisms as index permutations. Min RMSD over these is symmetry-
    corrected RMSD (what RDKit GetBestRMS computes): a molecule with equivalent atoms (symmetric
    substituents, ring pseudo-symmetry) must not be penalised for relabelling them. Falls back to
    the identity permutation if RDKit is unavailable or the graph can't be built. The heavy atoms of
    AddHs(MolFromSmiles) keep indices 0..n_heavy-1, so these indices align with the heavy-atom
    coordinate arrays used below."""
    try:
        from rdkit import Chem
    except ImportError:
        return [np.arange(n_heavy)]
    m = Chem.MolFromSmiles(smiles)
    if m is None or m.GetNumAtoms() != n_heavy:
        return [np.arange(n_heavy)]
    matches = m.GetSubstructMatches(m, uniquify=False, maxMatches=4096)
    return [np.array(t) for t in matches] or [np.arange(n_heavy)]


def rmsd_matrix_sym(A, B, perms):
    """Symmetry-corrected Kabsch RMSD: the min over automorphism permutations of B's atoms. With the
    identity permutation only, this is exactly `rmsd_matrix`."""
    out = None
    for p in perms:
        d = rmsd_matrix(A, B[:, p, :])
        out = d if out is None else np.minimum(out, d)
    return out


def recovery(D, taus):
    """Fraction of the reference set (columns of D) with a candidate (rows) within tau — the
    interpretable coverage convention. D[i,j] = distance(candidate_i, reference_j)."""
    nn = D.min(axis=0)
    return [float(np.mean(nn <= t)) for t in taus]


def tfd_full_matrix(smiles, n_atoms, coords_list):
    """All-pairs Torsion-Fingerprint-Deviation matrix over the given conformers (each n_atoms*3 in
    AddHs order). TFD (Schulz-Gasch 2012) is bounded [0,1], ring-aware and symmetry-aware — the
    field convention for macrocycle ensembles, where Cartesian RMSD conflates ring pseudo-rotation
    with global shape. Returns None if RDKit/TFD is unavailable or the molecule has no torsions
    (rigid)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import TorsionFingerprints
    except ImportError:
        return None
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if m is None or m.GetNumAtoms() != n_atoms:
        return None
    for c in coords_list:
        conf = Chem.Conformer(n_atoms)
        for i in range(n_atoms):
            conf.SetAtomPosition(i, (float(c[i][0]), float(c[i][1]), float(c[i][2])))
        m.AddConformer(conf, assignId=True)
    try:
        flat = TorsionFingerprints.GetTFDMatrix(m)
    except (RuntimeError, ValueError, IndexError):
        return None
    k = m.GetNumConformers()
    if k < 2 or len(flat) != k * (k - 1) // 2:
        return None
    M = np.zeros((k, k))
    idx = 0
    for i in range(k):
        for j in range(i):
            M[i, j] = M[j, i] = flat[idx]
            idx += 1
    return M


def bootstrap_ci(values, stat=np.median, n=3000, seed=20210185):
    """95% percentile-bootstrap CI of `stat` over molecules. Fixed seed → reproducible."""
    v = np.asarray(values, float)
    if len(v) == 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    boots = [stat(rng.choice(v, size=len(v), replace=True)) for _ in range(n)]
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def rmsd_matrix(A, B):
    """Kabsch-aligned RMSD between every conformer of A and every conformer of B."""
    Ac = A - A.mean(1, keepdims=True)
    Bc = B - B.mean(1, keepdims=True)
    n = A.shape[1]
    gA = (Ac ** 2).sum((1, 2))
    gB = (Bc ** 2).sum((1, 2))
    H = np.einsum("ikm,jkn->ijmn", Ac, Bc)
    U, S, Vt = np.linalg.svd(H)
    ds = np.sign(np.linalg.det(np.einsum("ijab,ijbc->ijac", U, Vt)))
    E = S[..., 0] + S[..., 1] + ds * S[..., 2]
    return np.sqrt(np.maximum((gA[:, None] + gB[None, :] - 2 * E) / n, 0.0))


# ---------------------------------------------------------------- stage 1

def stage1():
    """Perception: bounds matrix, elements, bonds, torsions, against RDKit's own values."""
    ref = json.load(open(PERCEPTION))
    smi = smiles_by_id()
    worst = 0.0
    failures = []
    for mid, r in sorted(ref.items()):
        o = spec_for(smi[mid])
        n = o["n_atoms"]
        bad = []
        if n != r["n_atoms"]:
            bad.append("n_atoms")
        if o["formal_charge"] != r["formal_charge"]:
            bad.append("formal_charge")
        if o["atomic_numbers"] != r["atomic_numbers"]:
            bad.append("elements")
        if sorted(sorted(b) for b in o["bonds"]) != [list(b) for b in r["bonds"]]:
            bad.append("bonds")
        d = max(abs(o["bounds"][i * n + j] - r["bounds"][i][j]) for i in range(n) for j in range(n))
        worst = max(worst, d)
        if d > BOUNDS_TOL:
            bad.append(f"bounds({d:.5f})")
        # RDKit's Python API reports only SMARTS-matched library torsions; the bridge additionally
        # carries basic-knowledge planarity terms. The relation is containment, not equality.
        def k(atoms, v):
            return (tuple(sorted(atoms)), tuple(round(float(x), 4) for x in v))
        ours = {k(t["atoms"], t["v"]) for t in o["exp_torsions"]}
        theirs = {k(t[:4], t[4:]) for t in r["torsions"]}
        if theirs - ours:
            bad.append(f"{len(theirs - ours)} RDKit torsions absent")
        if bad:
            failures.append((mid, bad))
    print(f"  stage1 perception : {len(ref)} molecules, worst bounds delta {worst:.6f} A")
    for mid, bad in failures[:8]:
        print(f"    FAIL {mid}: {bad}")
    return not failures


# ---------------------------------------------------------------- stage 2

def stage2():
    """Recipe parameters and pinned atoms.

    Recipe parameters come from her `count_exo_rotatable` extracted by AST and executed, which needs
    no embedding and so covers the whole corpus. Pinned atoms are only reachable by tracing her file,
    so that check is narrow by necessity.
    """
    smi = smiles_by_id()
    ok = True

    exo_ref = {}
    for p in sorted(FIX.glob(EXO_GLOB)):
        exo_ref.update(json.load(open(p)))
    if exo_ref:
        n_exo, bad_exo, exo_branches = 0, [], {}
        for name, r in exo_ref.items():
            mid = mol_id(name)
            if mid not in smi:
                continue
            n_exo += 1
            k = (r["core_seeds"], r["sidechain_confs"])
            exo_branches[k] = exo_branches.get(k, 0) + 1
            o = spec_for(smi[mid])
            if (o["core_seeds"], o["sidechain_confs"]) != k:
                bad_exo.append(mid)
        br = ", ".join(f"{a}x{b}: {v}" for (a, b), v in sorted(exo_branches.items()))
        print(f"  stage2 recipe     : {n_exo} molecules, {len(bad_exo)} mismatches [{br}]")
        for mid in bad_exo[:5]:
            print(f"    FAIL recipe {mid}")
        ok &= not bad_exo
    ref = {}
    sources = sorted(FIX.glob(STAGE2_GLOB))
    for p in sources:
        ref.update(json.load(open(p)))
    n = 0
    branches = {}
    bad_counts, bad_pins, varying = [], [], []
    for name, r in sorted(ref.items()):
        mid = mol_id(name)
        if mid not in smi:
            continue
        n += 1
        o = spec_for(smi[mid])
        branches[(r["rdkit_confs_1"], r["sidechain_confs"])] = \
            branches.get((r["rdkit_confs_1"], r["sidechain_confs"]), 0) + 1
        if (o["core_seeds"], o["sidechain_confs"]) != (r["rdkit_confs_1"], r["sidechain_confs"]):
            bad_counts.append(mid)
        if sorted(o["pin_atoms"]) != sorted(r["pins"]):
            bad_pins.append(mid)
        # her pin set is expected to be identical across a molecule's core seeds
        if r.get("pin_sets_distinct", 1) != 1:
            varying.append((mid, r["pin_sets_distinct"]))
    br = ", ".join(f"{a}x{b}: {k}" for (a, b), k in sorted(branches.items()))
    print(f"  stage2 pin atoms  : {n} molecules, {len(bad_pins)} mismatches "
          f"[traced branches — {br}]")
    for mid in bad_counts[:5]:
        print(f"    FAIL counts {mid} (traced)")
    for mid in bad_pins[:5]:
        print(f"    FAIL pins {mid}")
    for mid, k in varying[:5]:
        print(f"    note {mid}: oracle used {k} distinct pin sets across seeds")
    return ok and n > 0 and not bad_counts and not bad_pins


# ---------------------------------------------------------------- stage 3

def _wb(D, seeds, side):
    """within-seed (mean pairwise inside each sidechain block) and between-seed (mean pairwise among
    the block-leading core conformers) spread, from a full distance matrix D."""
    within = np.mean([
        D[j * side:(j + 1) * side, j * side:(j + 1) * side][np.triu_indices(side, 1)].mean()
        for j in range(seeds)
    ])
    cores = [j * side for j in range(seeds)]
    between = D[np.ix_(cores, cores)][np.triu_indices(seeds, 1)].mean()
    return float(within), float(between)


# Provisional acceptance band (STEP 4 will replace these with her own independent-seed run-to-run
# spread; until then they are the historical ±15% placeholders, and the gate is applied to the
# bootstrap CI, not a point estimate).
NN_MAX, DIV_MIN = 1.15, 0.85
TAUS = [0.5, 1.0, 1.5, 2.0]  # Angstrom, for the recovery curve

# --- Criterion-2 gate: the ONE definition, imported by criterion2_gate.py / fairness_control.py so the
#     verdict can never diverge by entry point. `load_band` is the ONE band loader with a single fallback.
PROD_SEED = 210185  # her production embed base (core seeds PROD_SEED+j, sidechain seed PROD_SEED)
REC_TOL = 0.9       # native must recover her confs to >= REC_TOL x her own half's recovery at each tau


def load_band():
    """The ratified acceptance band (`derived_band.json`) if present, else the single provisional
    fallback. One source for every gate entry point (check.stage3, criterion2_gate, fairness_control),
    so a missing band file yields the SAME threshold everywhere instead of a per-caller constant."""
    bp = FIX / "derived_band.json"
    if bp.exists():
        return json.load(open(bp))
    return {"NN_MAX": NN_MAX, "DIV_MIN": DIV_MIN, "source": "provisional (no derived_band.json)"}


def gate_failed(row, div_min):
    """Criterion-2 per-molecule verdict (the ratified relation — no under-coverage + no collapse): FAIL
    iff native under-covers her conformers (recovery below REC_TOL x her own half at any tau where her
    recovery is meaningful) or collapses diversity (within/between below the band floor `div_min`)."""
    collapse = row["within"] < div_min or row["between"] < div_min
    under = any(o < REC_TOL * h for h, o in zip(row["rec_her"], row["rec_our"]) if h >= 0.05)
    return collapse or under


def score(spec_txt, conformers):
    """3D Stage-A bounds-violation energy per conformer, via `bb-embed <spec> score < coords`. Writes the
    spec to a private tempdir (never a fixed repo path), so concurrent callers can't clobber each other."""
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        sp = pathlib.Path(td) / "spec.json"
        sp.write_text(spec_txt)
        payload = json.dumps({"conformers": [list(map(float, c)) for c in conformers]})
        out = subprocess.run([str(BB_EMBED), str(sp), "score"], input=payload,
                             stdout=subprocess.PIPE, universal_newlines=True).stdout
    return np.array(json.loads(out)["energies"], float)


def _ratio_row(smi_str, spec, Hf, Of, side, her_nconf):
    """One molecule's spread-normalised ratios. Hf = the reference (base-1) ensemble, Of = the
    CHALLENGER ensemble being tested for equivalence (native's embed in stage3; her independent
    second seed base in `band`). Both are flat coord arrays (conf, n_atoms, 3) in AddHs order.
    Returns (row, None) on success or (None, partial_tuple) if the ensembles can't be block-aligned.

    The reference is split into two independent seed-block halves HA|HB. Coverage of HB is measured
    two ways: by her own other half HA (the between-independent-seed baseline) and by the challenger.
    The ratio challenger/baseline is 1.0 when the challenger reproduces her ensemble exactly as well
    as an independent run of her own pipeline does — which is the target-B claim."""
    hi = heavy(spec["atomic_numbers"])
    na = spec["n_atoms"]
    seeds = spec["core_seeds"]
    # Either ensemble may be missing whole seed-blocks (a core-embed failure drops that seed's block),
    # but each surviving block is still complete (`side` confs). Compare on the ACTUAL block counts
    # (Kh reference, Ko challenger) rather than discarding the molecule; only a non-block-aligned
    # count or <2 blocks (can't split the reference into halves) is unusable.
    Kh, Ko = len(Hf) // side, len(Of) // side
    if len(Hf) % side or Kh < 2 or Ko * side < (Kh // 2) * side:
        return None, (mol_id(spec.get("name", "")) or "?", her_nconf, seeds * side)
    H, O = Hf[:, hi, :], Of[:, hi, :]
    half = (Kh // 2) * side
    perms = automorphisms(smi_str, len(hi))

    # symmetry-corrected RMSD matrices (rows = first arg, cols = second)
    D_HH = rmsd_matrix_sym(H, H, perms)
    D_OO = rmsd_matrix_sym(O, O, perms)
    D_OH = rmsd_matrix_sym(O, H, perms)
    her_nn = float(np.median(D_HH[:half, half:2 * half].min(axis=1)))   # HA vs HB
    our_nn = float(np.median(D_OH[:half, half:2 * half].min(axis=1)))   # O[:half] vs HB
    hw, hb = _wb(D_HH, Kh, side)
    ow, ob = _wb(D_OO, Ko, side)
    rec_her = recovery(D_HH[:half, half:2 * half], TAUS)
    rec_our = recovery(D_OH[:half, half:2 * half], TAUS)

    row = dict(her_nn=her_nn, nn=our_nn / her_nn if her_nn else float("nan"),
               within=ow / hw if hw else float("nan"), between=ob / hb if hb else float("nan"),
               rec_her=rec_her, rec_our=rec_our, tfd_nn=None, tfd_within=None, tfd_between=None)

    # TFD (ring-aware, symmetry-aware); one combined matrix over reference + challenger
    M = tfd_full_matrix(smi_str, na, list(Hf) + list(Of))
    if M is not None:
        n = len(Hf)
        T_HH, T_OO, T_OH = M[:n, :n], M[n:, n:], M[n:, :n]
        t_her = float(np.median(T_HH[:half, half:2 * half].min(axis=1)))
        t_our = float(np.median(T_OH[:half, half:2 * half].min(axis=1)))
        thw, thb = _wb(T_HH, Kh, side)
        tow, tob = _wb(T_OO, Ko, side)
        row.update(tfd_nn=t_our / t_her if t_her else float("nan"),
                   tfd_within=tow / thw if thw else float("nan"),
                   tfd_between=tob / thb if thb else float("nan"))
    return row, None


def stage3():
    """Ensembles vs the ORACLE's own seed-to-seed spread. Distances are symmetry-corrected RMSD
    (min over graph automorphisms) and TFD; coverage is reported both as a spread-normalised
    nearest-neighbour ratio and as an absolute recovery curve against her own recovery; diversity is
    guarded within- and between-seed; aggregates carry a bootstrap 95% CI and a per-molecule pass
    count, so a good average cannot hide a molecule that collapses."""
    # The acceptance band + gate criterion come from the single source (`load_band`/`gate_failed`), so
    # this and the standalone gate scripts can't diverge (see the module-level definitions above).
    band = load_band()
    nn_max, div_min, band_src = band["NN_MAX"], band["DIV_MIN"], band.get("source", "?")
    smi = smiles_by_id()
    rows = []
    partial = []
    for mid, smiles, spec, her, Hf, _f in iter_ensembles(ENSEMBLES, smi):
        side = spec["sidechain_confs"]
        o = embed(spec, PROD_SEED)
        Of = np.array(o["conformers"]).reshape(-1, o["n_atoms"], 3)
        row, part = _ratio_row(smiles, spec, Hf, Of, side, her["n_conformers"])
        if row is None:
            partial.append((mid, part[1], part[2]))
            continue
        row["mid"] = mid
        rows.append(row)

    print(f"  stage3 ensembles  : {len(rows)} molecules  (symmetry-corrected RMSD + TFD)")
    for mid, got, want in partial:
        print(f"    skipped {mid}: oracle produced {got} of {want} conformers (seed blocks unknown)")

    # RATIFIED relation (criterion 2): a molecule FAILS iff native UNDER-covers her conformers or
    # COLLAPSES diversity — the benchmark's two named failure modes. Native being MORE diverse (the
    # symmetric NN-band "fail") is NOT a failure mode. The verdict is `gate_failed` (module-level, the
    # ONE definition, also used by the standalone gate scripts); see [[conformance-contract]].
    def failed(r):
        return gate_failed(r, div_min)

    hdr = f"    {'mol':14s} {'nn':>5s} {'within':>6s} {'betwn':>6s} | {'tfd_nn':>6s} {'t_wi':>5s} {'t_bt':>5s}"
    print(hdr)
    for r in rows:
        flag = " FAIL" if failed(r) else ""
        t = lambda k: f"{r[k]:.2f}" if r[k] is not None and not math.isnan(r[k]) else "  -"
        print(f"    {r['mid']:14s} {r['nn']:5.2f} {r['within']:6.2f} {r['between']:6.2f} | "
              f"{t('tfd_nn'):>6s} {t('tfd_within'):>5s} {t('tfd_between'):>5s}{flag}")

    def agg(key):
        vals = [r[key] for r in rows if r[key] is not None and not math.isnan(r[key])]
        lo, hi = bootstrap_ci(vals)
        return float(np.median(vals)), lo, hi, len(vals)

    npass = sum(not failed(r) for r in rows)
    print(f"    -- aggregates (median [95% CI], n) --")
    for key, lbl in [("nn", "RMSD NN"), ("within", "RMSD within"), ("between", "RMSD between"),
                     ("tfd_nn", "TFD NN"), ("tfd_within", "TFD within"), ("tfd_between", "TFD between")]:
        m, lo, hi, n = agg(key)
        print(f"      {lbl:12s} {m:5.2f}x  [{lo:.2f}, {hi:.2f}]   (n={n})")
    rh = np.mean([r["rec_her"] for r in rows], axis=0)
    ro = np.mean([r["rec_our"] for r in rows], axis=0)
    print(f"    -- recovery of her HB (fraction within tau) --")
    print(f"      tau (A):   " + "  ".join(f"{t:5.1f}" for t in TAUS))
    print(f"      her  half: " + "  ".join(f"{x:5.2f}" for x in rh))
    print(f"      ours     : " + "  ".join(f"{x:5.2f}" for x in ro))
    print(f"    per-molecule pass (RATIFIED: no under-coverage + no collapse; DIV_MIN={div_min:.2f}, "
          f"REC_TOL={REC_TOL}): {npass}/{len(rows)}")
    for r in rows:
        if failed(r):
            print(f"      FAIL {r['mid']}: within={r['within']:.2f} between={r['between']:.2f} "
                  f"rec_her={['%.2f'%x for x in r['rec_her']]} rec_our={['%.2f'%x for x in r['rec_our']]}")

    # RATIFIED gate: no under-coverage (mean recovery of her >= REC_TOL x her own half at every tau) AND
    # no collapse (diversity CI not below her run-to-run floor) AND >=95% of molecules pass per-molecule.
    # The NN band is reported above for information but is NOT a gate criterion (it penalizes valid
    # enrichment, not a criterion-2 failure mode).
    rec_ok = all(o >= REC_TOL * h for h, o in zip(rh, ro) if h >= 0.05)
    _, wi_lo, _, _ = agg("within")
    _, bt_lo, _, _ = agg("between")
    return rec_ok and wi_lo >= div_min and bt_lo >= div_min and npass >= 0.95 * len(rows)


def band():
    """Derive the acceptance band EMPIRICALLY, from her own between-independent-seed-base spread.

    stage3 normalises native's coverage of her ensemble by her OWN half-to-half coverage, so a ratio
    of 1.0 is perfect. But how far from 1.0 is still 'as good as her own pipeline'? That tolerance is
    not a guess: run her verbatim recipe at a SECOND, fully independent seed base (SEED_BASE != 210185,
    generated by gen_reference.py into BB_ENSEMBLES2) and compute the SAME ratios with base-2 standing
    in for native. The resulting distribution IS her run-to-run variability; its upper (for coverage)
    and lower (for diversity) bootstrap-CI edges are the band stage3 must fall inside. This replaces
    the historical ±15% placeholder with a measured quantity (CONFORMANCE.md target B, STEP 4)."""
    ens2 = os.environ.get("BB_ENSEMBLES2")
    if not ens2:
        print("  band              : skipped (set BB_ENSEMBLES2=<dir> to a second-seed-base run)")
        return True
    D2 = pathlib.Path(ens2) if os.path.isabs(ens2) else FIX / ens2
    smi = smiles_by_id()
    rows, partial = [], []
    for mid, smiles, spec, her, Hf, f in iter_ensembles(ENSEMBLES, smi):
        g = D2 / f.name
        if not g.exists():
            continue
        her2 = json.load(open(g))
        side = spec["sidechain_confs"]
        Of = np.array(her2["conformers"]).reshape(-1, her2["n_atoms"], 3)
        row, part = _ratio_row(smiles, spec, Hf, Of, side, her["n_conformers"])
        if row is None:
            partial.append((mid, part[1], part[2]))
            continue
        row["mid"] = mid
        rows.append(row)
    if not rows:
        print(f"  band              : no molecules present in both {ENSEMBLES.name} and {D2.name}")
        return True

    def agg(key):
        vals = [r[key] for r in rows if r[key] is not None and not math.isnan(r[key])]
        lo, hi = bootstrap_ci(vals)
        return float(np.median(vals)), lo, hi, len(vals)

    print(f"  band (base-2)     : {len(rows)} molecules  (her {ENSEMBLES.name} vs her {D2.name})")
    print(f"    -- her between-seed-base ratios (median [95% CI], n) --")
    derived = {}
    for key in ("nn", "within", "between", "tfd_nn", "tfd_within", "tfd_between"):
        m, lo, hi, n = agg(key)
        derived[key] = dict(median=m, lo=lo, hi=hi, n=n)
        print(f"      {key:12s} {m:5.2f}x  [{lo:.2f}, {hi:.2f}]   (n={n})")
    # The band: coverage ratios (nn) may be as high as her own base-to-base upper CI; diversity ratios
    # (within/between) as low as her own lower CI. Native must fall inside. A symmetric guard is used
    # for nn (some slack below 1 too) but the binding edge is the upper one.
    out = dict(NN_MAX=round(derived["nn"]["hi"], 3),
               DIV_MIN=round(min(derived["within"]["lo"], derived["between"]["lo"]), 3),
               TFD_NN_MAX=round(derived["tfd_nn"]["hi"], 3) if not math.isnan(derived["tfd_nn"]["hi"]) else None,
               source=f"{ENSEMBLES.name} vs {D2.name}", n=len(rows), ratios=derived)
    bp = FIX / "derived_band.json"
    json.dump(out, open(bp, "w"), indent=1)
    print(f"    -> derived band: NN_MAX={out['NN_MAX']} DIV_MIN={out['DIV_MIN']} "
          f"TFD_NN_MAX={out['TFD_NN_MAX']}  (wrote {bp.name}; stage3 will consume it)")
    return True


# ---------------------------------------------------------------- stage 4

def read_mol2_coords(path):
    """Flat xyz from a mol2 ATOM block, in file order."""
    coords, inblock = [], False
    for line in open(path):
        if line.startswith("@<TRIPOS>ATOM"):
            inblock = True
            continue
        if line.startswith("@<TRIPOS>") and inblock:
            break
        if inblock and line.strip():
            p = line.split()
            coords += [float(p[2]), float(p[3]), float(p[4])]
    return coords


def stage4():
    """Solvation: her geometry in, her `.solv` out.

    Her solvation conformer comes from an unseeded embed (`build_ligands.py:112`), so an end-to-end
    comparison would be distributional. Feeding bb-solv the geometry she actually gave AMSOL makes
    the whole chain deterministic — AMSOL on fixed coordinates is reproducible — so any difference
    is ours: Z-matrix construction, input generation, invocation, or output parsing.

    Needs BB_ORACLE_RUN (a completed run) and the AMSOL settings bb-solv reads from the environment.
    """
    run = os.environ.get("BB_ORACLE_RUN")
    if not run or not BB_SOLVATE.exists():
        print("  stage4 solvation  : skipped (set BB_ORACLE_RUN; needs bb-solvate)")
        return True
    run = pathlib.Path(run)
    # a completed run may be given as the working directory or as its parent
    if not (run / "solv").is_dir() and (run / "out" / "solv").is_dir():
        run = run / "out"
    smi = smiles_by_id()
    identical, differing, skipped = [], [], []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        for solv_path in sorted(run.glob("solv/*/output.solv")):
            name = solv_path.parent.name.replace(".mol2", "")
            mid = mol_id(name)
            mol2 = run / "3d" / f"{name}.mol2"
            if mid not in smi or not mol2.exists():
                skipped.append(name)
                continue
            spec = spec_for(smi[mid])
            coords = read_mol2_coords(mol2)
            if len(coords) // 3 != spec["n_atoms"]:
                skipped.append(f"{name} (atom count)")
                continue
            sp = tmp / f"{name}.spec.json"
            cf = tmp / f"{name}.confs.json"
            sp.write_text(json.dumps(spec))
            cf.write_text(json.dumps({"n_atoms": len(coords) // 3, "conformers": [coords]}))
            out = tmp / name
            r = subprocess.run([str(BB_SOLVATE), str(sp), str(cf), str(out), "0"],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
            ours = out / "output.solv"
            if r.returncode != 0 or not ours.exists():
                differing.append((name, (r.stderr or "bb-solvate failed").strip()[:70]))
                continue
            if ours.read_text() == solv_path.read_text():
                identical.append(name)
            else:
                differing.append((name, "output.solv differs"))
    print(f"  stage4 solvation  : {len(identical)} byte-identical, {len(differing)} differing, "
          f"{len(skipped)} skipped")
    for name, why in differing[:5]:
        print(f"    FAIL {name}: {why}")
    return bool(identical) and not differing


# ---------------------------------------------------------------- rings

def rings():
    """Ring perception in pure Rust against RDKit's symmetrized SSSR.

    Two things are checked, and the second is the requirement. Exact ring-set equality is reported
    because it is the sharpest signal, but equal-size alternates in bridged bicyclics are chosen by
    RDKit's traversal rather than by a rule, so exact equality is not the contract. What must hold is
    that every quantity the pipeline reads from the ring set agrees: the ring-atom union behind
    `count_exo_rotatable`, largest-ring selection for the pins, and ring-size membership, which is
    what `r{9-}` in the patched macrocycle torsions keys on.
    """
    ref_path = FIX / os.environ.get("BB_RINGS", "rings_seeds5000.json")
    if not ref_path.exists() or not BB_PERCEIVE.exists():
        print("  rings perception  : skipped (needs rings fixture and bb-perceive)")
        return True
    ref = json.load(open(ref_path))
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_PERCEIVE), str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "rings" in d:
            ours[d["name"]] = [sorted(r) for r in d["rings"]]

    def derived(rs):
        union, sizes = set(), {}
        for r in rs:
            for a in r:
                union.add(a)
                sizes.setdefault(a, set()).add(len(r))
        largest = set(max(rs, key=len)) if rs else set()
        r9 = {a for r in rs if len(r) >= 9 for a in r}
        return union, largest, r9, sizes

    n = exact = 0
    bad_union, bad_largest, bad_r9 = [], [], []
    for name, r in ref.items():
        o = ours.get(name)
        if o is None:
            continue
        n += 1
        R = sorted(sorted(x) for x in r["rings"])
        O = sorted(o)
        if R == O:
            exact += 1
        u1, l1, r91, _ = derived(R)
        u2, l2, r92, _ = derived(O)
        if u1 != u2:
            bad_union.append(name)
        if l1 != l2:
            bad_largest.append(name)
        if r91 != r92:
            bad_r9.append(name)
    print(f"  rings perception  : {n} molecules, {exact} exact ({100.0 * exact / max(n, 1):.2f}%)")
    print(f"    consumed quantities — ring-atom union {len(bad_union)}, "
          f"largest ring {len(bad_largest)}, r{{9-}} membership {len(bad_r9)} mismatches")
    for name in (bad_union + bad_largest + bad_r9)[:5]:
        print(f"    FAIL {name}")
    return n > 0 and not bad_union and not bad_largest and not bad_r9


def hydrogens():
    """Valence and hydrogen counts in pure Rust against RDKit's, per atom.

    Compared as sequences rather than totals: a molecule whose hydrogen count is right in aggregate
    but wrong per atom yields a different AddHs ordering, and every later index would be misaligned.
    """
    ref_path = FIX / os.environ.get("BB_HCOUNT", "hcount_seeds5000.json")
    if not ref_path.exists() or not BB_PERCEIVE.exists():
        print("  hydrogens         : skipped (needs hcount fixture and bb-perceive)")
        return True
    ref = json.load(open(ref_path))
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_PERCEIVE), str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "total_num_hs" in d:
            ours[d["name"]] = d
    n = atoms = 0
    bad = {"explicit_valence": [], "implicit_valence": [], "total_num_hs": []}
    for name, r in ref.items():
        o = ours.get(name)
        if o is None:
            continue
        n += 1
        atoms += r["n_atoms"]
        for field in bad:
            if o[field] != r[field]:
                bad[field].append(name)
    print(f"  hydrogens         : {n} molecules, {atoms} atoms")
    for field, names in bad.items():
        print(f"    {field:18s} {len(names)} mismatches")
    for field, names in bad.items():
        for name in names[:3]:
            print(f"    FAIL {field} {name}")
    return n > 0 and not any(bad.values())


def addhs():
    """Explicit hydrogens: atom order and bond order after AddHs, against RDKit's.

    Both are compared as sequences. Atom order is what every per-atom array downstream means; bond
    order is read by ring perception, and a bond list that matches as a set but not as a sequence
    silently changes which of two equal-size rings is found.
    """
    ref_path = FIX / os.environ.get("BB_ADDHS", "addhs_seeds5000.json")
    if not ref_path.exists() or not BB_PERCEIVE.exists():
        print("  addhs             : skipped (needs addhs fixture and bb-perceive)")
        return True
    ref = json.load(open(ref_path))
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_PERCEIVE), str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "addhs_z" in d:
            ours[d["name"]] = d
    n = atoms = 0
    bad_atoms, bad_seq, bad_set = [], [], []
    for name, r in ref.items():
        o = ours.get(name)
        if o is None:
            continue
        n += 1
        atoms += r["n_atoms"]
        if o["addhs_z"] != r["z"]:
            bad_atoms.append(name)
        ob = [tuple(x) for x in o["addhs_bonds"]]
        rb = [tuple(x) for x in r["bonds"]]
        if ob != rb:
            bad_seq.append(name)
        if sorted(tuple(sorted(x)) for x in ob) != sorted(tuple(sorted(x)) for x in rb):
            bad_set.append(name)
    print(f"  addhs             : {n} molecules, {atoms} atoms")
    print(f"    atom sequence      {len(bad_atoms)} mismatches")
    print(f"    bond sequence      {len(bad_seq)} mismatches "
          f"(as a set: {len(bad_set)})")
    for name in (bad_atoms + bad_seq)[:4]:
        print(f"    FAIL {name}")
    return n > 0 and not bad_atoms and not bad_seq


def aromaticity():
    """The SMILES's written aromaticity against RDKit's perceived aromaticity.

    The front-end does not perceive aromaticity; it reads what the SMILES states. That is only valid
    while the input is written in aromatic form — Kekule-form input (`C1=CC=CC=C1`) writes no
    aromatic atom and RDKit perceives six, so the two would diverge completely. This gate is what
    makes that an assumption under test rather than a silent dependency on how the corpus was
    generated. If it ever fails, RDKit's aromaticity model has to be implemented.
    """
    ref_path = FIX / os.environ.get("BB_KEKULIZE", "kekulize_seeds5000.json")
    if not ref_path.exists() or not BB_PERCEIVE.exists():
        print("  aromaticity       : skipped (needs kekulize fixture and bb-perceive)")
        return True
    ref = json.load(open(ref_path))
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_PERCEIVE), str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "arom_written" in d:
            ours[d["name"]] = d
    n = atoms = arom = 0
    bad = []
    for name, r in ref.items():
        o = ours.get(name)
        if o is None:
            continue
        n += 1
        atoms += r["n_atoms"]
        arom += sum(r["aromatic_atoms"])
        if o["arom_written"] != r["aromatic_atoms"]:
            bad.append(name)
    print(f"  aromaticity       : {n} molecules, {atoms} atoms ({arom} aromatic), "
          f"{len(bad)} mismatches")
    for name in bad[:4]:
        print(f"    FAIL {name} — written aromaticity differs from perceived")
    return n > 0 and not bad


def hybrid():
    """Conjugation and hybridization in Rust against RDKit's.

    Hybridization reaches beyond its obvious uses: UFF atom typing keys on it, so it feeds bond
    lengths and angles in `setTopolBounds`, and inversion terms apply only at SP2 centres. It is
    computed over heavy atoms; RDKit sanitizes before `AddHs`, so appended hydrogens are left
    unspecified.
    """
    ref_path = FIX / os.environ.get("BB_HYBRID", "hybrid_seeds5000.json")
    if not ref_path.exists() or not BB_PERCEIVE.exists():
        print("  hybridization     : skipped (needs hybrid fixture and bb-perceive)")
        return True
    ref = json.load(open(ref_path))
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_PERCEIVE), str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "hybridization" in d:
            ours[d["name"]] = d
    n = atoms = conj = 0
    bad_h, bad_c = [], []
    for name, r in ref.items():
        o = ours.get(name)
        if o is None:
            continue
        n += 1
        atoms += r["n_atoms"]
        conj += sum(r["conjugated"])
        if o["hybridization"] != r["hybridization"]:
            bad_h.append(name)
        if o["conjugated"] != r["conjugated"]:
            bad_c.append(name)
    print(f"  hybridization     : {n} molecules, {atoms} atoms, {len(bad_h)} mismatches")
    print(f"  conjugation       : {conj} conjugated bonds, {len(bad_c)} mismatches")
    for name in (bad_h + bad_c)[:4]:
        print(f"    FAIL {name}")
    return n > 0 and not bad_h and not bad_c


def angles_check():
    """Bond angles in Rust against RDKit's collectBondsAndAngles.

    Compared as a sequence: the angle list is generated by iterating bond pairs, so its order is
    inherited from bond order and a set comparison would not notice the difference.
    """
    ref_path = FIX / os.environ.get("BB_ANGLES", "angles_seeds600.json")
    if not ref_path.exists() or not BB_PERCEIVE.exists():
        print("  angles            : skipped (needs angles fixture and bb-perceive)")
        return True
    ref = json.load(open(ref_path))
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_PERCEIVE), str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "angles" in d:
            ours[d["name"]] = d
    n = tot = 0
    bad_seq, bad_set = [], []
    for name, r in ref.items():
        o = ours.get(name)
        if o is None:
            continue
        n += 1
        tot += len(r["angles"])
        ra = [tuple(x) for x in r["angles"]]
        oa = [tuple(x) for x in o["angles"]]
        if oa != ra:
            bad_seq.append(name)
        if sorted(oa) != sorted(ra):
            bad_set.append(name)
    print(f"  angles            : {n} molecules, {tot} angles, {len(bad_seq)} sequence "
          f"mismatches (as a set: {len(bad_set)})")
    for name in bad_seq[:4]:
        print(f"    FAIL {name}")
    return n > 0 and not bad_seq


def recipe():
    """Recipe parameters and the pinned core, computed in Rust from our own ring perception.

    The same values stage 2 checks, but produced by `bb-perceive` rather than the RDKit bridge —
    `count_exo_rotatable` against her AST-extracted function over the whole corpus, and the largest
    ring against the traced coordMap.
    """
    if not BB_PERCEIVE.exists():
        print("  recipe (Rust)     : skipped (needs bb-perceive)")
        return True
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_PERCEIVE), str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "exo" in d:
            ours[d["name"]] = d

    exo_ref = {}
    for p in sorted(FIX.glob(EXO_GLOB)):
        exo_ref.update(json.load(open(p)))
    n = 0
    bad_exo, bad_recipe = [], []
    for name, r in exo_ref.items():
        o = ours.get(name)
        if o is None:
            continue
        n += 1
        if o["exo"] != r["exo"]:
            bad_exo.append(name)
        if (o["core_seeds"], o["sidechain_confs"]) != (r["core_seeds"], r["sidechain_confs"]):
            bad_recipe.append(name)
    print(f"  recipe (Rust)     : {n} molecules, {len(bad_exo)} exo, "
          f"{len(bad_recipe)} recipe mismatches")

    pins = {}
    for p in sorted(FIX.glob(STAGE2_GLOB)):
        pins.update(json.load(open(p)))
    m, bad_core = 0, []
    for name, r in pins.items():
        o = ours.get(mol_id(name))
        if o is None:
            continue
        m += 1
        if sorted(o["core_atoms"]) != sorted(r["core_atoms"]):
            bad_core.append(mol_id(name))
    print(f"  pinned core (Rust): {m} molecules, {len(bad_core)} mismatches")
    for name in (bad_exo + bad_recipe + bad_core)[:4]:
        print(f"    FAIL {name}")
    return n > 0 and not bad_exo and not bad_recipe and not bad_core


def uff_labels():
    """UFF atom typing in pure Rust against RDKit's getAtomLabel.

    The label decides which UFF parameter row set12Bounds reads (r1 -> bond rest length), so it is
    the correctness foundation of the whole bounds matrix. The Rust port reproduces the full
    getAtomLabel + addAtomChargeFlags source (every element branch); this gate is exact per-atom
    equality over the corpus, not the source of truth for what to implement.
    """
    ref_path = FIX / "uff_labels_seeds600.jsonl"
    if not ref_path.exists() or not BB_BOUNDS_DUMP.exists():
        print("  uff atom typing   : skipped (needs uff_labels fixture and bb-uff-labels)")
        return True
    ref = {}
    for line in open(ref_path):
        d = json.loads(line)
        if "labels" in d:
            ref[d["name"]] = d["labels"]
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_BOUNDS_DUMP), "labels", str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "labels" in d:
            ours[d["name"]] = d["labels"]

    n = atoms = mism = 0
    bad = []
    for name, R in ref.items():
        O = ours.get(name)
        if O is None:
            continue
        n += 1
        if len(R) != len(O):
            bad.append((name, f"len {len(R)}!={len(O)}"))
            continue
        for x, y in zip(R, O):
            atoms += 1
            if x != y:
                mism += 1
                if len(bad) < 8:
                    bad.append((name, f"{y}!={x}"))
    print(f"  uff atom typing   : {n} molecules, {atoms} atoms, {mism} label mismatches")
    for name, why in bad[:5]:
        print(f"    FAIL {name}: {why}")
    return n > 0 and mism == 0 and not any("len" in w for _, w in bad)


def uff_bondlen():
    """UFF bond rest lengths in pure Rust against RDKit's calcBondRestLength.

    These are set12Bounds' accumData.bondLengths — the 1-2 distance every bound in the matrix is
    anchored to. The Rust port reproduces both branches of set12: calcBondRestLength when both atoms
    have UFF params, and RDKit's (rvdw1+rvdw2)/2 fallback otherwise. The gate is bit-exact equality
    over the corpus. Over this corpus every atom has params, so the fallback branch is not exercised
    here — it is gated separately by the uff module's dummy-atom unit test.
    """
    ref_path = FIX / "uff_bondlen_seeds600.jsonl"
    if not ref_path.exists() or not BB_BOUNDS_DUMP.exists():
        print("  uff bond lengths  : skipped (needs uff_bondlen fixture and bb-uff-bondlen)")
        return True
    ref = {}
    for line in open(ref_path):
        d = json.loads(line)
        if "bondlen" in d:
            ref[d["name"]] = d["bondlen"]
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_BOUNDS_DUMP), "bondlen", str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "bondlen" in d:
            ours[d["name"]] = d["bondlen"]

    n = bonds = 0
    worst = 0.0
    bad = []
    for name, R in ref.items():
        O = ours.get(name)
        if O is None or len(O) != len(R):
            bad.append((name, "missing/len"))
            continue
        n += 1
        for x, y in zip(R, O):
            bonds += 1
            d = abs(x - y)
            worst = max(worst, d)
            if d > 1e-9 and len(bad) < 8:
                bad.append((name, f"{y} vs {x}"))
    print(f"  uff bond lengths  : {n} molecules, {bonds} bonds, worst |delta| {worst:.2e}")
    for name, why in bad[:5]:
        print(f"    FAIL {name}: {why}")
    return n > 0 and worst == 0.0 and not bad


def topo_dist():
    """Topological distance matrix in pure Rust against RDKit's getDistanceMat.

    setLowerBoundVDW reads this for the 1-5 (dist==4) and 1-6 (dist==5) van der Waals scaling, and
    set14 uses it. The gate is exact per-pair equality over the corpus. This is more thorough than
    the transitive coverage the bounds gate gives, which only exercises the dist==4/dist==5
    distinctions.
    """
    import gzip
    ref_path = FIX / "topo_dist_seeds600.jsonl.gz"
    if not ref_path.exists() or not BB_BOUNDS_DUMP.exists():
        print("  topo distance     : skipped (needs topo_dist fixture and bb-bounds-dump)")
        return True
    ref = {}
    for line in gzip.open(ref_path, "rt"):
        d = json.loads(line)
        if "dist" in d:
            ref[d["name"]] = d["dist"]
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_BOUNDS_DUMP), "dist", str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if "dist" in d:
            ours[d["name"]] = d["dist"]

    n = pairs = mism = 0
    bad = []
    for name, R in ref.items():
        O = ours.get(name)
        if O is None or len(O) != len(R):
            bad.append((name, "missing/len"))
            continue
        n += 1
        for x, y in zip(R, O):
            pairs += 1
            if x != y:
                mism += 1
                if len(bad) < 8:
                    bad.append((name, f"{y} vs {x}"))
    print(f"  topo distance     : {n} molecules, {pairs} pairs, {mism} mismatches")
    for name, why in bad[:5]:
        print(f"    FAIL {name}: {why}")
    return n > 0 and mism == 0 and not bad


def _bounds_stage_gate(mode, fixture, tol, label, key=None):
    """Shared driver for a staged bounds-matrix gate (bounds12, bounds13, ...).

    `mode` is the bb-bounds-dump subcommand; `key` is the JSON field both sides use (defaults to
    `mode`, but the full matrix reuses the raw_bounds fixture whose field is "bounds").
    """
    import gzip
    key = key or mode
    ref_path = FIX / fixture
    if not ref_path.exists() or not BB_BOUNDS_DUMP.exists():
        print(f"  {label}: skipped (needs {fixture} and bb-bounds-dump)")
        return True
    ref = {}
    opener = gzip.open if str(ref_path).endswith(".gz") else open
    for line in opener(ref_path, "rt"):
        d = json.loads(line)
        if key in d:
            ref[d["name"]] = d[key]
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_BOUNDS_DUMP), mode, str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {}
    for line in out.splitlines():
        d = json.loads(line)
        if key in d:
            ours[d["name"]] = d[key]
    n = cells = 0
    worst = 0.0
    bad = []
    for name, R in ref.items():
        O = ours.get(name)
        if O is None or len(O) != len(R):
            bad.append((name, "missing/len"))
            continue
        n += 1
        for x, y in zip(R, O):
            cells += 1
            d = abs(x - y)
            worst = max(worst, d)
            if d > tol and len(bad) < 8:
                bad.append((name, f"{y} vs {x}"))
    print(f"  {label}: {n} molecules, {cells} cells, worst |delta| {worst:.2e}")
    for name, why in bad[:5]:
        print(f"    FAIL {name}: {why}")
    return n > 0 and worst <= tol and not bad


def bounds12():
    """set12 + setLowerBoundVDW vs RDKit staged setTopolBounds — bit-exact."""
    return _bounds_stage_gate("bounds12", "bounds12_seeds600.jsonl.gz", 0.0, "bounds12(set12+vdw)")


def bounds13():
    """set12 + set13 + VDW vs RDKit — to float precision (set13 uses cos/sqrt, so not bit-exact)."""
    return _bounds_stage_gate("bounds13", "bounds13_seeds600.jsonl.gz", 1e-9, "bounds13(+set13) ")


def bounds14():
    """set12 + set13 + set14 + VDW vs RDKit — 1-4 torsion layer incl. macrocycle + double-bond
    stereo. Float precision (compute14Dist* uses cos/sin/sqrt)."""
    return _bounds_stage_gate("bounds14", "bounds14_seeds600.jsonl.gz", 1e-6, "bounds14(+set14) ")


def bounds15():
    """The full setTopolBounds (set12+13+14+15+VDW) vs RDKit's raw (pre-smoothing) bounds — the
    complete pure-Rust bounds matrix. Float precision (compute15Dists uses cos/sin/sqrt/acos)."""
    return _bounds_stage_gate("bounds15", "raw_bounds_seeds600.jsonl.gz", 1e-6,
                              "bounds15(full)  ", key="bounds")


def chiral_tags():
    """Tetrahedral chirality in pure Rust against RDKit's getChiralTag.

    Reproduces AdjustAtomChiralityFlags (@/@@ + GetBondOrdering perturbation + tag inversion). Exact
    per-atom equality; only the tetrahedral tags (0/1/2) are compared, non-tetrahedral types skipped.
    """
    import gzip
    ref_path = FIX / "chiral_tags_seeds600.jsonl.gz"
    if not ref_path.exists() or not BB_BOUNDS_DUMP.exists():
        print("  chiral tags       : skipped (needs chiral_tags fixture and bb-bounds-dump)")
        return True
    ref = {}
    for line in gzip.open(ref_path, "rt"):
        d = json.loads(line)
        if "tags" in d:
            ref[d["name"]] = d["tags"]
    corpus = next((p for p in CORPUS if "5000" in p.name), CORPUS[0] if CORPUS else None)
    out = subprocess.run([str(BB_BOUNDS_DUMP), "chiraltags", str(corpus)], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    ours = {json.loads(l)["name"]: json.loads(l).get("tags") for l in out.splitlines()}
    n = atoms = mism = 0
    bad = []
    for name, R in ref.items():
        O = ours.get(name)
        if O is None or len(O) != len(R):
            bad.append((name, "missing/len"))
            continue
        n += 1
        for r, o in zip(R, O):
            if r not in (0, 1, 2):  # skip non-tetrahedral chiral types
                continue
            atoms += 1
            if r != o:
                mism += 1
                if len(bad) < 8:
                    bad.append((name, f"{o}!={r}"))
    print(f"  chiral tags       : {n} molecules, {atoms} atoms, {mism} tag mismatches")
    for name, why in bad[:5]:
        print(f"    FAIL {name}: {why}")
    return n > 0 and mism == 0 and not any("len" in w for _, w in bad)


STAGES = {"rings": rings, "hydrogens": hydrogens, "addhs": addhs, "aromaticity": aromaticity,
          "recipe": recipe, "angles": angles_check, "hybrid": hybrid, "uff_labels": uff_labels,
          "uff_bondlen": uff_bondlen, "topo_dist": topo_dist, "chiral_tags": chiral_tags,
          "bounds12": bounds12, "bounds13": bounds13, "bounds14": bounds14, "bounds15": bounds15,
          "stage1": stage1, "stage2": stage2, "stage3": stage3, "band": band, "stage4": stage4}

if __name__ == "__main__":
    want = [a for a in sys.argv[1:] if a in STAGES] or list(STAGES)
    results = {name: STAGES[name]() for name in want}
    print("\n  " + "  ".join(f"{k}={'PASS' if v else 'FAIL'}" for k, v in results.items()))
    sys.exit(0 if all(results.values()) else 1)
