#!/usr/bin/env python3
# Does native's EXTRA diversity (the conformers beyond her distribution — the source of the ~2x TFD
# richness) consist of VALID conformers or JUNK? This is the forcing check for criterion 2: if native's
# excess confs are as bounds-satisfying as its in-distribution confs, the enrichment is valid coverage
# (Accept is forced); if the excess is high-strain, it's a real quality regression (Hold).
#
# For each molecule: native confs, each labelled by (a) nearest-neighbour RMSD to HER ensemble
# (heavy-atom, Kabsch) — "excess" = far from her; and (b) its 3D Stage-A bounds-violation energy
# (bb-embed score, the shared certified spec both pipelines target). Then compare the bounds energy of
# the EXCESS confs (top-quartile distance-to-her) vs the IN-DISTRIBUTION confs (bottom-quartile).
#
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/parity/excess_validity.py [N]
import json
import subprocess
import sys
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check  # single source: smiles_by_id, mol_id, native_spec_for, embed, score, rmsd_matrix, PROD_SEED

ENS = check.REPO / "validation/fixtures/divya/ensembles"


def main():
    n_mol = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    smi = check.smiles_by_id()
    rows = []
    for f in sorted(ENS.glob("*.confs.json"))[:n_mol]:
        her = json.load(open(f))
        mid = check.mol_id(her["name"])
        if mid not in smi:
            continue
        spec = check.native_spec_for(smi[mid])
        na = spec["n_atoms"]
        if her["n_atoms"] != na:
            continue
        hi = [i for i, z in enumerate(spec["atomic_numbers"]) if z > 1]
        nat = check.embed(spec, check.PROD_SEED)
        Nf = np.array(nat["conformers"]).reshape(-1, na, 3)
        Hf = np.array(her["conformers"]).reshape(-1, na, 3)
        # bounds-violation energy per native conf (shared certified spec)
        be = check.score(json.dumps(spec), nat["conformers"])
        # nearest-her RMSD per native conf (heavy atoms)
        D = check.rmsd_matrix(Nf[:, hi, :], Hf[:, hi, :])
        nn = D.min(axis=1)
        # split native confs by distance-to-her: excess = top quartile, core = bottom quartile
        q1, q3 = np.percentile(nn, 25), np.percentile(nn, 75)
        core = be[nn <= q1]; excess = be[nn >= q3]
        if len(core) and len(excess):
            rows.append((mid, float(nn.mean()), float(core.mean()), float(excess.mean()),
                         float(np.median(be)), float(be.max())))
            print(f"  {mid:14s} nn-to-her {rows[-1][1]:.2f}A | bounds-E core {rows[-1][2]:.2f} "
                  f"excess {rows[-1][3]:.2f}  (ratio {rows[-1][3]/max(rows[-1][2],1e-6):.2f})")
            sys.stdout.flush()
    if not rows:
        print("no molecules"); return 1
    core = np.array([r[2] for r in rows]); excess = np.array([r[3] for r in rows])
    ratio = excess / np.maximum(core, 1e-6)
    print(f"\n=== {len(rows)} molecules: are native's EXCESS (far-from-her) confs valid or junk? ===")
    print(f"  bounds-violation energy (shared certified spec):")
    print(f"    IN-DISTRIBUTION (near her, bottom-quartile dist): median-of-means {np.median(core):.3f}")
    print(f"    EXCESS (far from her, top-quartile dist):          median-of-means {np.median(excess):.3f}")
    print(f"    excess/core ratio: median {np.median(ratio):.2f}  mean {ratio.mean():.2f}  (>>1 = junk, ~1 = valid)")
    print(f"  => native's extra diversity is {'VALID coverage (Accept forced)' if np.median(ratio) < 1.5 else 'ELEVATED-strain (revisit)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
