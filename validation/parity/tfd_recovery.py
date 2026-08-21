#!/usr/bin/env python3
# The decisive "no under-coverage on the ring-sensitive metric" measure for criterion 2. The stage3
# TFD-NN ratio (native->her) conflates native's extra diversity with under-coverage. This measures the
# clean direction: TFD RECOVERY (her->native) = fraction of HER conformers that have a NATIVE conformer
# within a TFD threshold, compared against her OWN other half's recovery of her (the reference). If
# native recovers her ring-torsion conformers at parity, "no under-coverage" holds on TFD and native's
# TFD-NN 1.42 is benign extra-diversity (Accept forced). If native under-recovers, it's a real shortfall.
#
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/parity/tfd_recovery.py [N]
import sys
import pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check  # reuse tfd_full_matrix, native_spec_for, embed, smiles_by_id, mol_id, ENSEMBLES

TAUS = [0.05, 0.10, 0.15, 0.20]  # TFD is in [0,1]; these are typical ring-torsion recovery radii


def main():
    n_mol = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    smi = check.smiles_by_id()
    her_rec, our_rec = [], []
    n = 0
    for mid, smiles, spec, her, Hf, f in check.iter_ensembles(check.ENSEMBLES, smi):
        if n >= n_mol:
            break
        na = spec["n_atoms"]
        side = spec["sidechain_confs"]
        Kh = len(Hf) // side
        if her["n_atoms"] != na or Kh < 2:
            continue
        o = check.embed(spec, check.PROD_SEED)
        Of = np.array(o["conformers"]).reshape(-1, o["n_atoms"], 3)
        half = (Kh // 2) * side
        M = check.tfd_full_matrix(smiles, na, list(Hf) + list(Of))
        if M is None:
            continue
        nH = len(Hf)
        T_HH = M[:nH, :nH]
        T_OH = M[nH:, :nH]  # rows = native (candidates), cols = her
        # recover HER second half (cols half:2half) by: her OWN first half (HA) vs NATIVE
        her_rec.append(check.recovery(T_HH[:half, half:2 * half], TAUS))  # HA covers HB
        our_rec.append(check.recovery(T_OH[:, half:2 * half], TAUS))       # native covers HB
        n += 1
        print(f"  {mid:14s} her-half {['%.2f'%x for x in her_rec[-1]]}  native {['%.2f'%x for x in our_rec[-1]]}")
        sys.stdout.flush()
    if not her_rec:
        print("no molecules"); return 1
    hr = np.mean(her_rec, axis=0); orr = np.mean(our_rec, axis=0)
    print(f"\n=== TFD recovery of HER conformers ({n} molecules) — no-under-coverage test ===")
    print(f"  tau (TFD):        " + "  ".join("%5.2f" % t for t in TAUS))
    print(f"  by HER OWN half:  " + "  ".join("%5.2f" % x for x in hr))
    print(f"  by NATIVE:        " + "  ".join("%5.2f" % x for x in orr))
    ratio = orr / np.maximum(hr, 1e-9)
    print(f"  native/her ratio: " + "  ".join("%5.2f" % r for r in ratio))
    ok = np.all(ratio >= 0.9)  # native recovers her within 10% of her own half at every radius
    print(f"\n  => native {'RECOVERS her ring-torsion confs at parity (no under-coverage on TFD -> Accept forced)' if ok else 'UNDER-recovers her ring-torsion confs (real shortfall -> Hold)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
