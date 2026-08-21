#!/usr/bin/env python3
# Bounds-strain comparison (goal criterion 2, strain axis) on the RIGHT ruler. MMFF (strain_compare.py)
# is a poor absolute ruler here -- neither pipeline optimizes against it and it scores her own
# production conformers at ~57 kcal/mol strain. The decision-relevant geometric ruler is the 3D
# Stage-A BOUNDS-VIOLATION energy: how well a conformer satisfies the DistGeom bounds matrix that BOTH
# pipelines embed to (certified byte-identical, [[conformance-contract]]). native's Stage-A minimizer
# optimizes exactly this. Scoring her ensemble and native's on the SAME spec bounds asks: are native's
# (more dispersed) conformers as geometrically valid, by the shared spec, as hers?
#
# RDKit-FREE: bb-spec-native + bb-embed only -> runs host-side.  validation/parity/bounds_strain.py
import json
import subprocess
import sys
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check  # single source: smiles_by_id, mol_id, native_spec_for, embed, score, PROD_SEED

ENS = check.REPO / "validation/fixtures/divya/ensembles"


def main():
    smi = check.smiles_by_id()
    rows = []
    for f in sorted(ENS.glob("*.confs.json")):
        her = json.load(open(f))
        name = her["name"]
        mid = check.mol_id(name)
        if mid not in smi:
            continue
        spec = check.native_spec_for(smi[mid])
        na = spec["n_atoms"]
        if her["n_atoms"] != na:
            continue
        nat = check.embed(spec, check.PROD_SEED)
        Hf = np.array(her["conformers"]).reshape(-1, na * 3)
        Nf = np.array(nat["conformers"]).reshape(-1, na * 3)
        hE = check.score(json.dumps(spec), Hf)
        nE = check.score(json.dumps(spec), Nf)
        if len(hE) == 0 or len(nE) == 0:
            continue
        rows.append(dict(mid=mid,
                         h_med=float(np.median(hE)), n_med=float(np.median(nE)),
                         h_p90=float(np.percentile(hE, 90)), n_p90=float(np.percentile(nE, 90)),
                         h_max=float(hE.max()), n_max=float(nE.max()),
                         frac_above=float(np.mean(nE > hE.max()))))
        r = rows[-1]
        print(f"  {mid:14s} bounds-E med h/n {r['h_med']:5.2f}/{r['n_med']:5.2f}  "
              f"p90 {r['h_p90']:5.2f}/{r['n_p90']:5.2f}  max {r['h_max']:6.2f}/{r['n_max']:6.2f}  "
              f"native>her_max {r['frac_above']*100:4.1f}%")
        sys.stdout.flush()

    if not rows:
        print("no molecules scored")
        return 1

    def med(k):
        return float(np.median([r[k] for r in rows]))

    print("\n=== aggregate over %d molecules (3D Stage-A bounds-violation energy, shared certified bounds) ===" % len(rows))
    print(f"  median bounds-E:  her {med('h_med'):.2f}   native {med('n_med'):.2f}   (ratio {med('n_med')/med('h_med'):.2f})")
    print(f"  p90    bounds-E:  her {med('h_p90'):.2f}   native {med('n_p90'):.2f}")
    print(f"  max    bounds-E:  her {med('h_max'):.2f}   native {med('n_max'):.2f}")
    fa = float(np.median([r["frac_above"] for r in rows]))
    fam = float(np.mean([r["frac_above"] for r in rows]))
    print(f"  JUNK: native confs above HER max bounds-E:  median {fa*100:.1f}%  mean {fam*100:.1f}%")
    worst = sorted(rows, key=lambda r: -r["frac_above"])[:6]
    print("  most-junk:", ", ".join(f"{r['mid']}({r['frac_above']*100:.0f}%)" for r in worst))
    json.dump(rows, open(pathlib.Path(__file__).parent / "bounds_strain.json", "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
