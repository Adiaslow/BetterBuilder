#!/usr/bin/env python3
# Fairness control for the ratified stage3 gate. The gate FAILS native (29/96 per-molecule). Is that a
# real shortfall, or is the gate so strict that even HER OWN independent run can't clear it? This runs
# the IDENTICAL gate (no under-coverage: recovery >= REC_TOL x her own half at each tau; no collapse:
# diversity >= DIV_MIN) with her SECOND seed base (base2) standing in for native. If base2 passes ~all,
# the gate is fair and native's failure is real; if base2 also fails heavily, the gate is too strict.
#
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/parity/fairness_control.py
import json
import sys
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check

REPO = HERE.parent.parent  # validation/parity -> validation -> repo root
B1 = REPO / "validation/fixtures/divya/ensembles"
B2 = REPO / "validation/fixtures/divya/ensembles_base2"
# gate criterion + band loader come from check — the ONE source (docs/2026-08-21-code-audit.md C3).


def main():
    smi = check.smiles_by_id()
    div_min = check.load_band()["DIV_MIN"]
    npass = n = 0
    fails = []
    for mid, smiles, spec, her, Hf, f in check.iter_ensembles(B1, smi):
        g = B2 / f.name
        if not g.exists():
            continue
        her2 = json.load(open(g))
        side = spec["sidechain_confs"]
        Of = np.array(her2["conformers"]).reshape(-1, her2["n_atoms"], 3)
        row, _ = check._ratio_row(smiles, spec, Hf, Of, side, her["n_conformers"])
        if row is None:
            continue
        n += 1
        if check.gate_failed(row, div_min):
            fails.append((mid, row))
        else:
            npass += 1
    print(f"=== FAIRNESS CONTROL: her OWN base2 run vs base1, SAME ratified gate (DIV_MIN={div_min:.2f}, REC_TOL={check.REC_TOL}) ===")
    print(f"  base2 (her independent run) per-molecule pass: {npass}/{n}")
    print(f"  (native scored 29/96 on the same gate)")
    if n:
        print(f"  => the gate is {'FAIR (her own run passes; native genuinely under-covers)' if npass >= 0.9 * n else 'TOO STRICT (even her own independent run fails it) -> native failure is within her run-to-run'}")
    for mid, r in fails[:6]:
        print(f"    base2 FAIL {mid}: within={r['within']:.2f} between={r['between']:.2f} "
              f"rec_her={['%.2f'%x for x in r['rec_her']]} rec_our={['%.2f'%x for x in r['rec_our']]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
