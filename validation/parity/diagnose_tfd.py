#!/usr/bin/env python3
# TFD-gap diagnosis (dim-2 investigation): stage3 shows native reproduces her ensemble in ABSOLUTE
# recovery but with a TFD-NN ratio ~1.44x vs her ~1.0 between-seed floor. This decomposes that gap
# PER TORSION to answer three things, so the acceptance-criterion decision rests on mechanism:
#   1. Is it a few molecules or uniform?  (per-molecule tfd_nn distribution)
#   2. SHIFT vs DISPERSION: does native land at a different rotamer (basin difference) or explore a
#      WIDER angle range (over-dispersion)?  per torsion: circular-mean shift + circular-std ratio.
#   3. RING vs EXOCYCLIC: is the divergence in the macrocycle backbone (ring torsions, hard) or in
#      rotatable exocyclic single bonds (likely a fixable minimizer/filter difference)?
#
# Runs in her container (needs RDKit TorsionFingerprints):
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/parity/diagnose_tfd.py
import json
import math
import os
import subprocess
import sys
import pathlib
import numpy as np
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints, rdMolTransforms

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check  # single source: smiles_by_id, mol_id, native_spec_for, embed, PROD_SEED
REPO = check.REPO
ENS = REPO / "validation/fixtures/divya/ensembles"
# CHALLENGER: what plays the role compared against her base-1 ensemble.
#   "native" (default) -> native embed_recipe;  a dir path -> her ensemble from that dir (e.g. base-2).
# The base-2 run is the HER-VS-HER baseline: same shift/disp metrics with her own second seed base,
# so native's numbers can be read against her own between-seed floor rather than against 0.
CHALLENGER = os.environ.get("DIAG_CHALLENGER", "native")


def native_ensemble(smiles):
    d = check.embed(check.native_spec_for(smiles), check.PROD_SEED)
    return np.array(d["conformers"]).reshape(-1, d["n_atoms"], 3), d["n_atoms"]


def set_conf(mol, coords):
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))
    cid = mol.AddConformer(conf, assignId=True)
    return cid


def circ_stats(angles_deg):
    """circular mean (deg) and circular std (deg) of a set of dihedral angles."""
    a = np.deg2rad(np.asarray(angles_deg, float))
    C, S = np.cos(a).mean(), np.sin(a).mean()
    R = math.hypot(C, S)
    mean = math.degrees(math.atan2(S, C))
    std = math.degrees(math.sqrt(-2.0 * math.log(R))) if R > 1e-12 else 180.0
    return mean, std


def circ_dist(a, b):
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


def torsion_angles(mol, quads, cids):
    """for each torsion (first quad), the dihedral (deg) across the given conformer ids."""
    out = []
    for cid in cids:
        conf = mol.GetConformer(cid)
        out.append([rdMolTransforms.GetDihedralDeg(conf, *q) for q in quads])
    return np.array(out)  # (n_conf, n_tors)


def main():
    smi = check.smiles_by_id()
    # per-molecule tfd_nn from the stage3 log, if available, to correlate
    tfd_nn = {}
    logp = sys.argv[2] if len(sys.argv) > 2 else os.environ.get("BB_STAGE3_LOG", "")
    if os.path.exists(logp):
        for line in open(logp):
            p = line.split()
            if len(p) >= 8 and p[0].startswith("CSLB") and p[4] == "|":
                try:
                    tfd_nn[p[0]] = float(p[5])
                except ValueError:
                    pass

    rows = []
    for f in sorted(ENS.glob("*.confs.json")):
        her = json.load(open(f))
        name = her["name"]
        mid = check.mol_id(name)
        if mid not in smi:
            continue
        m = Chem.AddHs(Chem.MolFromSmiles(smi[mid]))
        na = m.GetNumAtoms()
        Hf = np.array(her["conformers"]).reshape(-1, her["n_atoms"], 3)
        if her["n_atoms"] != na:
            continue
        try:
            if CHALLENGER == "native":
                Nf, nna = native_ensemble(smi[mid])
            else:
                g = pathlib.Path(CHALLENGER)
                g = g if g.is_absolute() else REPO / g
                d2 = json.load(open(g / f.name))
                Nf, nna = np.array(d2["conformers"]).reshape(-1, d2["n_atoms"], 3), d2["n_atoms"]
        except Exception:
            continue
        if nna != na:
            continue
        # torsion definitions (non-ring + ring), each entry (list_of_quads, weight); is_ring flag
        try:
            nonring, ring = TorsionFingerprints.CalculateTorsionLists(m)
        except Exception:
            continue
        tors = [(t[0][0], False) for t in nonring if t[0]] + [(t[0][0], True) for t in ring if t[0]]
        if not tors:
            continue
        quads = [q for q, _ in tors]
        is_ring = np.array([r for _, r in tors])

        # build conformers on a scratch mol for her and native separately
        mh = Chem.Mol(m); mh.RemoveAllConformers()
        hc = [set_conf(mh, c) for c in Hf]
        mn = Chem.Mol(m); mn.RemoveAllConformers()
        nc = [set_conf(mn, c) for c in Nf]
        HA = torsion_angles(mh, quads, hc)   # (nH, T)
        NA = torsion_angles(mn, quads, nc)   # (nN, T)

        shifts, disp = [], []
        for t in range(len(quads)):
            hm, hs = circ_stats(HA[:, t])
            nm, ns = circ_stats(NA[:, t])
            shifts.append(circ_dist(nm, hm))               # basin shift (deg)
            disp.append(ns / hs if hs > 1e-6 else float("nan"))  # dispersion ratio (native/her)
        shifts = np.array(shifts); disp = np.array(disp)
        rows.append(dict(mid=mid, T=len(quads), nring=int(is_ring.sum()),
                         shift_ring=float(np.nanmean(shifts[is_ring])) if is_ring.any() else float("nan"),
                         shift_exo=float(np.nanmean(shifts[~is_ring])) if (~is_ring).any() else float("nan"),
                         disp_ring=float(np.nanmean(disp[is_ring])) if is_ring.any() else float("nan"),
                         disp_exo=float(np.nanmean(disp[~is_ring])) if (~is_ring).any() else float("nan"),
                         tfd_nn=tfd_nn.get(mid, float("nan"))))
        print(f"  {mid:14s} T={len(quads):2d} ring={int(is_ring.sum()):2d} | "
              f"shift ring/exo {rows[-1]['shift_ring']:5.1f}/{rows[-1]['shift_exo']:5.1f}deg  "
              f"disp ring/exo {rows[-1]['disp_ring']:4.2f}/{rows[-1]['disp_exo']:4.2f}  "
              f"tfd_nn={rows[-1]['tfd_nn']:.2f}")
        sys.stdout.flush()

    if not rows:
        print("no molecules diagnosed")
        return 1

    def mean(key):
        v = [r[key] for r in rows if not math.isnan(r[key])]
        return float(np.mean(v)) if v else float("nan")

    print("\n=== aggregate over %d molecules (challenger=%s) ===" % (len(rows), CHALLENGER))
    print(f"  basin SHIFT (native mean vs her mean, deg):   ring {mean('shift_ring'):5.1f}   exo {mean('shift_exo'):5.1f}")
    print(f"  DISPERSION ratio (native std / her std):       ring {mean('disp_ring'):4.2f}    exo {mean('disp_exo'):4.2f}")
    print("  (shift>>0 => different rotamer basin; disp>>1 => native over-explores; <1 => under-explores)")

    # correlate per-molecule dispersion/shift with tfd_nn to see what drives the gap
    good = [r for r in rows if not math.isnan(r["tfd_nn"])]
    if good:
        import numpy as _np
        y = _np.array([r["tfd_nn"] for r in good])
        for key in ("disp_exo", "disp_ring", "shift_exo", "shift_ring"):
            x = _np.array([r[key] for r in good])
            ok = ~_np.isnan(x)
            if ok.sum() > 5 and _np.std(x[ok]) > 1e-9:
                c = float(_np.corrcoef(x[ok], y[ok])[0, 1])
                print(f"  corr(tfd_nn, {key:9s}) = {c:+.2f}")
        hi = sorted(good, key=lambda r: -r["tfd_nn"])[:8]
        print("  worst-TFD molecules:", ", ".join(f"{r['mid']}({r['tfd_nn']:.1f})" for r in hi))
    tag = "native" if CHALLENGER == "native" else "base2"
    json.dump(rows, open(HERE / ("tfd_diag_%s.json" % tag), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
