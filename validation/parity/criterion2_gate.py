#!/usr/bin/env python3
# THE ratified criterion-2 gate (distributional reproduction), CALIBRATED to her own run-to-run.
#
# The absolute "per-molecule pass >= 95%" bar is unachievable even by HER OWN pipeline at any strictness
# that resolves a real difference (her base2-vs-base1 passes only ~18% of the strict recovery+diversity
# per-molecule test — 100-conf recovery fractions fluctuate more than the tolerance between her own
# reseeds). So "statistically indistinguishable within her own run-to-run variability" is operationalised
# correctly as: native's per-molecule pass rate on the IDENTICAL strict test is >= her own independent
# run's (base2). If native is at least as self-consistent with her ensemble as she is with herself, it is
# within her run-to-run band. Failable: if a native-embed regression made it under-cover MORE than her own
# reseed does, native_pass < base2_pass and the gate goes RED.
#
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/parity/criterion2_gate.py [N]
import json
import sys
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check

REPO = HERE.parent.parent
B1 = REPO / "validation/fixtures/divya/ensembles"
B2 = REPO / "validation/fixtures/divya/ensembles_base2"
# gate criterion (REC_TOL, load_band, gate_failed) and PROD_SEED come from check — the ONE source, so
# this gate can't diverge from check.stage3 or fairness_control (see docs/2026-08-21-code-audit.md C3).


def main():
    n_mol = int(sys.argv[1]) if len(sys.argv) > 1 else 96
    smi = check.smiles_by_id()
    div_min = check.load_band()["DIV_MIN"]
    nat_pass = nat_n = b2_pass = b2_n = seen = 0
    for mid, smiles, spec, her, Hf, f in check.iter_ensembles(B1, smi):
        if seen >= n_mol:
            break
        side = spec["sidechain_confs"]
        if her["n_atoms"] != spec["n_atoms"]:
            continue
        seen += 1
        # native challenger
        o = check.embed(spec, check.PROD_SEED)
        Of = np.array(o["conformers"]).reshape(-1, o["n_atoms"], 3)
        row, _ = check._ratio_row(smiles, spec, Hf, Of, side, her["n_conformers"])
        if row is not None:
            nat_n += 1
            nat_pass += 0 if check.gate_failed(row, div_min) else 1
        # her base2 challenger (the calibration)
        g = B2 / f.name
        if g.exists():
            her2 = json.load(open(g))
            O2 = np.array(her2["conformers"]).reshape(-1, her2["n_atoms"], 3)
            row2, _ = check._ratio_row(smiles, spec, Hf, O2, side, her["n_conformers"])
            if row2 is not None:
                b2_n += 1
                b2_pass += 0 if check.gate_failed(row2, div_min) else 1
        sys.stdout.flush()
    nat_rate = nat_pass / max(nat_n, 1)
    b2_rate = b2_pass / max(b2_n, 1)
    print("=== CRITERION 2 (calibrated to her run-to-run) ===")
    print(f"  strict per-molecule test (no under-coverage REC_TOL={check.REC_TOL} + no collapse DIV_MIN={div_min:.2f}):")
    print(f"    native  : {nat_pass}/{nat_n}  ({100*nat_rate:.0f}%)")
    print(f"    her base2: {b2_pass}/{b2_n}  ({100*b2_rate:.0f}%)   <- her own run-to-run reference")
    ok = nat_rate >= b2_rate
    print(f"  => native {'>=' if ok else '<'} her own run-to-run  ->  criterion 2 {'PASS (within her variability)' if ok else 'FAIL (worse than her own reseed)'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
