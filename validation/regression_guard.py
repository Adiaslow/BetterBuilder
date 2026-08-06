"""Automated fidelity regression guard for the BetterBuilder candidate (durable, repo-relative).

Runs the candidate (fixed seed → deterministic) on the pinned reference molecules, compares its
conformer ensembles to the patched-oracle reference ensembles (heavy-atom NN-RMSD), and checks
per-molecule + aggregate metrics against a pinned baseline with noise-derived gates.

  <rdkit-python> validation/regression_guard.py --update-baseline   # re-pin the baseline
  <rdkit-python> validation/regression_guard.py                     # check; exit 0=PASS, 1=REGRESSION

Must run under a Python with RDKit (e.g. the patched engine's `oracle_python`); it shells out to the
release `bb-spec`/`bb-embed` binaries. Reference ensembles live in validation/reference/, the molecule
list in validation/molecules.smi, so the guard is self-contained apart from the built binaries.

Gates (from the measured 3-seed run-to-run noise: median NN 1.26-1.28, <2A mean 97%, <1.5A mean
69-70%, >=90%-count 88-92): FAIL if aggregate median NN > 1.32, or <2A mean < 95%, or <1.5A mean
< 64%, or >=90%-count < 85, or >=3 molecules worse by > 0.40 A median NN, or any mol worse by > 0.70 A.
"""
import sys, os, glob, json, subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem
from rdkit.Geometry import Point3D

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VALID = os.path.join(REPO, "validation")
REFDIR = os.path.join(VALID, "reference")
MOLECULES = os.path.join(VALID, "molecules.smi")
BASELINE = os.path.join(VALID, "regression_baseline.json")
BBSPEC = os.path.join(REPO, "rust/target/release/bb-spec")
# BB_EMBED_BIN overrides the embed binary (e.g. the GPU pipeline bb-embed-gpu) for a like-for-like check.
BBEMBED = os.environ.get("BB_EMBED_BIN", os.path.join(REPO, "rust/target/release/bb-embed"))
SEED = "210185"

# deterministic + no thread oversubscription (python workers parallelize; each bb-embed serial)
os.environ["RAYON_NUM_THREADS"] = "1"
os.environ["BB_SEED"] = SEED
# Hermetic by default: ignore a stray BB_GRAD_TOL in the environment (pass --allow-env to keep it,
# e.g. to deliberately inject a regression when self-testing the guard).
if "--allow-env" not in sys.argv:
    os.environ.pop("BB_GRAD_TOL", None)

# mcNNNN -> SMILES (self-contained)
name2smi = {}
for line in open(MOLECULES):
    p = line.split()
    if len(p) >= 2:
        name2smi[p[1].split("_")[0]] = p[0]


def all_pairs_rmsd(A, B):
    Ac = A - A.mean(1, keepdims=True); Bc = B - B.mean(1, keepdims=True)
    N = A.shape[1]; gA = (Ac**2).sum((1, 2)); gB = (Bc**2).sum((1, 2))
    H = np.einsum("ikm,jkn->ijmn", Ac, Bc); U, S, Vt = np.linalg.svd(H)
    ds = np.sign(np.linalg.det(np.einsum("ijab,ijbc->ijac", U, Vt)))
    E = S[..., 0] + S[..., 1] + ds * S[..., 2]
    return np.sqrt(np.maximum((gA[:, None] + gB[None, :] - 2 * E) / N, 0.0))


def _heavy(m):
    cf = m.GetConformer(); idx = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    return np.array([list(cf.GetAtomPosition(i)) for i in idx])


def metrics_for(name, smi, d, sdf):
    """NN-RMSD metrics for one molecule given its embed output dict d (n_atoms, conformers)."""
    n = d["n_atoms"]
    base = Chem.AddHs(Chem.MolFromSmiles(smi))
    if base.GetNumAtoms() != n:
        return None
    heavy_idx = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() > 1]
    A = np.stack([np.array([[c[3*i], c[3*i+1], c[3*i+2]] for i in heavy_idx]) for c in d["conformers"]])
    oracle = [m for m in Chem.SDMolSupplier(sdf, removeHs=False, sanitize=True) if m]
    B = np.stack([_heavy(m) for m in oracle])
    if A.shape[1] != B.shape[1]:
        return None
    D = all_pairs_rmsd(A, B); nn = np.concatenate([D.min(1), D.min(0)])
    return dict(med=float(np.median(nn)), cov20=float(np.mean(nn <= 2.0)), cov15=float(np.mean(nn <= 1.5)))


def one(job):
    """CPU-candidate path: embed one molecule in its own process (correct for the multi-core CPU)."""
    name, sdf = job
    smi = name2smi.get(name)
    if not smi:
        return (name, None)
    try:
        spec = subprocess.run([BBSPEC, smi], capture_output=True, text=True, timeout=120)
        if spec.returncode != 0:
            return (name, None)
        sf = os.path.join(VALID, f".rg_{name}.spec.json"); open(sf, "w").write(spec.stdout)
        emb = subprocess.run([BBEMBED, sf, "recipe", os.environ["BB_SEED"]],
                             capture_output=True, text=True, timeout=300)
        os.remove(sf)
        return (name, metrics_for(name, smi, json.loads(emb.stdout), sdf))
    except Exception:
        return (name, None)


def measure():
    jobs = [(os.path.basename(f).split("_")[0], f) for f in sorted(glob.glob(REFDIR + "/*.sdf"))]
    res = {}
    with ProcessPoolExecutor(max_workers=min(8, os.cpu_count())) as ex:
        for name, r in ex.map(one, jobs):
            if r is not None:
                res[name] = r
    return res


def _spec_job(job):
    name, sdf = job
    smi = name2smi.get(name)
    if not smi:
        return None
    p = subprocess.run([BBSPEC, smi], capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        return None
    return (name, smi, sdf, json.loads(p.stdout))


def measure_batch():
    """GPU path: specs on CPU (parallel), then ALL molecules embedded in ONE bb-embed-gpu process (one
    GPU context/command stream — no cross-process contention), then RMSD on CPU (parallel)."""
    jobs = [(os.path.basename(f).split("_")[0], f) for f in sorted(glob.glob(REFDIR + "/*.sdf"))]
    with ProcessPoolExecutor(max_workers=min(8, os.cpu_count())) as ex:
        specs = [s for s in ex.map(_spec_job, jobs) if s is not None]
    bf = os.path.join(VALID, ".rg_batch.jsonl")
    with open(bf, "w") as f:
        for name, smi, sdf, spec in specs:
            f.write(json.dumps({"name": name, "spec": spec}) + "\n")
    # ONE GPU process for every molecule.
    emb = subprocess.run([BBEMBED, bf, "batch", os.environ["BB_SEED"]], capture_output=True, text=True)
    os.remove(bf)
    by_name = {}
    for line in emb.stdout.splitlines():
        if line.strip():
            d = json.loads(line); by_name[d["name"]] = d
    meta = {name: (smi, sdf) for name, smi, sdf, _ in specs}
    args = [(name, meta[name][0], by_name[name], meta[name][1]) for name in by_name if name in meta]
    res = {}
    with ProcessPoolExecutor(max_workers=min(8, os.cpu_count())) as ex:
        for name, r in zip([a[0] for a in args], ex.map(_metrics_star, args)):
            if r is not None:
                res[name] = r
    return res


def _metrics_star(a):
    try:
        return metrics_for(*a)
    except Exception:
        return None


def aggregate(res):
    meds = np.array([r["med"] for r in res.values()])
    c20 = np.array([r["cov20"] for r in res.values()])
    c15 = np.array([r["cov15"] for r in res.values()])
    return dict(median_nn=float(np.median(meds)), cov20_mean=float(c20.mean()),
                cov15_mean=float(c15.mean()), cov90_count=int((c20 >= 0.9).sum()), n=len(res))


def main():
    # BB_EMBED_BATCH=1 → single-process GPU batch (one device context, no contention). Default = per-molecule
    # multiprocess (correct for the multi-core CPU candidate).
    res = measure_batch() if os.environ.get("BB_EMBED_BATCH") == "1" else measure()
    if not res:
        print("ERROR: no molecules measured (check binaries / reference data)"); return 2
    agg = aggregate(res)
    if "--update-baseline" in sys.argv:
        json.dump({"per_mol_med": {k: v["med"] for k, v in res.items()}, "agg": agg},
                  open(BASELINE, "w"), indent=1)
        print(f"baseline pinned: {agg['n']} molecules, median NN {agg['median_nn']:.3f}, "
              f"<2A {agg['cov20_mean']*100:.0f}%, >=90% {agg['cov90_count']}")
        return 0

    base = json.load(open(BASELINE)); bpm = base["per_mol_med"]
    fails = []
    if agg["median_nn"] > 1.32: fails.append(f"aggregate median NN {agg['median_nn']:.3f} > 1.32")
    if agg["cov20_mean"] < 0.95: fails.append(f"<2A mean {agg['cov20_mean']*100:.0f}% < 95%")
    if agg["cov15_mean"] < 0.64: fails.append(f"<1.5A mean {agg['cov15_mean']*100:.0f}% < 64%")
    if agg["cov90_count"] < 85: fails.append(f">=90% count {agg['cov90_count']} < 85")
    worse = sorted(((k, res[k]["med"] - bpm[k]) for k in bpm if k in res and res[k]["med"] - bpm[k] > 0.40),
                   key=lambda x: -x[1])
    if len(worse) >= 3:
        fails.append(f"{len(worse)} mols worse by >0.40A: " + ", ".join(f"{k}(+{d:.2f})" for k, d in worse[:6]))
    catastrophic = [(k, d) for k, d in worse if d > 0.70]
    if catastrophic:
        fails.append("catastrophic single-mol: " + ", ".join(f"{k}(+{d:.2f})" for k, d in catastrophic))

    print(f"current: {agg['n']} mols, median NN {agg['median_nn']:.3f} (base {base['agg']['median_nn']:.3f}), "
          f"<2A {agg['cov20_mean']*100:.0f}% (base {base['agg']['cov20_mean']*100:.0f}%), "
          f">=90% {agg['cov90_count']} (base {base['agg']['cov90_count']})")
    if fails:
        print("REGRESSION:")
        for f in fails:
            print("  - " + f)
        return 1
    print("PASS - within noise band of baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
