#!/usr/bin/env python3
# Strain/energy-distribution comparison (goal criterion 2, the un-measured axis). The TFD diagnosis
# showed native samples the CORRECT basins but is ~2x more torsionally dispersed than her ensemble.
# This asks the decisive follow-up: is native's extra diversity LOW-strain (benign -- more valid
# low-energy basins) or HIGH-strain (junk -- native emitting bad conformers)?
#
# Ruler: MMFF94 -- a NEUTRAL third-party energy neither pipeline optimized against (her conformers are
# ETKDG/DistGeom-optimized; native's are Stage-A/B/C-optimized), so it scores both ensembles without
# favouring either. UFF fallback for molecules MMFF can't type. Per molecule we compare native's and
# her per-conformer energy distributions (absolute, and relative to each ensemble's own minimum), and
# the high-strain tail (fraction of native confs above HER ensemble's max energy = the junk indicator).
#
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/parity/strain_compare.py
import json
import math
import os
import subprocess
import sys
import pathlib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check  # single source: smiles_by_id, mol_id, native_spec_for, embed, PROD_SEED
REPO = check.REPO
ENS = REPO / "validation/fixtures/divya/ensembles"


def native_ensemble(smiles):
    d = check.embed(check.native_spec_for(smiles), check.PROD_SEED)
    return np.array(d["conformers"]).reshape(-1, d["n_atoms"], 3), d["n_atoms"]


def energies(smiles, flat):
    """MMFF94 (UFF fallback) energy of each conformer (kcal/mol). Returns (list, ff_name)."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    na = m.GetNumAtoms()
    if flat.shape[1] != na:
        return None, None
    m.RemoveAllConformers()
    for c in flat:
        conf = Chem.Conformer(na)
        for i in range(na):
            conf.SetAtomPosition(i, (float(c[i][0]), float(c[i][1]), float(c[i][2])))
        m.AddConformer(conf, assignId=True)
    ff_name = "MMFF"
    props = AllChem.MMFFGetMoleculeProperties(m)
    out = []
    for cid in range(m.GetNumConformers()):
        ff = None
        if props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(m, props, confId=cid)
        if ff is None:
            ff = AllChem.UFFGetMoleculeForceField(m, confId=cid)
            ff_name = "UFF"
        if ff is None:
            return None, None
        out.append(ff.CalcEnergy())
    return np.array(out, float), ff_name


def main():
    smi = check.smiles_by_id()
    rows = []
    for f in sorted(ENS.glob("*.confs.json")):
        her = json.load(open(f))
        name = her["name"]
        mid = check.mol_id(name)
        if mid not in smi:
            continue
        Hf = np.array(her["conformers"]).reshape(-1, her["n_atoms"], 3)
        try:
            Nf, nna = native_ensemble(smi[mid])
        except Exception:
            continue
        if nna != her["n_atoms"]:
            continue
        hE, ff1 = energies(smi[mid], Hf)
        nE, ff2 = energies(smi[mid], Nf)
        if hE is None or nE is None or ff1 != ff2:
            continue
        # absolute per-conformer energy distributions and the strain (relative to each ensemble's min)
        h_med, n_med = float(np.median(hE)), float(np.median(nE))
        h_p90, n_p90 = float(np.percentile(hE, 90)), float(np.percentile(nE, 90))
        h_str, n_str = hE - hE.min(), nE - nE.min()      # strain above each ensemble's own minimum
        frac_above = float(np.mean(nE > hE.max()))        # native confs worse than HER worst = junk
        rows.append(dict(mid=mid, ff=ff1, n=len(nE),
                         d_med=n_med - h_med, d_p90=n_p90 - h_p90,
                         h_strmed=float(np.median(h_str)), n_strmed=float(np.median(n_str)),
                         h_strp90=float(np.percentile(h_str, 90)), n_strp90=float(np.percentile(n_str, 90)),
                         frac_above=frac_above, n_minus_h_min=float(nE.min() - hE.min())))
        r = rows[-1]
        print(f"  {mid:14s} {ff1} dMed {r['d_med']:+6.1f} dP90 {r['d_p90']:+6.1f} | "
              f"strain(med h/n) {r['h_strmed']:5.1f}/{r['n_strmed']:5.1f}  "
              f"(p90 {r['h_strp90']:5.1f}/{r['n_strp90']:5.1f})  junk {r['frac_above']*100:4.1f}%")
        sys.stdout.flush()

    if not rows:
        print("no molecules scored")
        return 1

    def med(key):
        return float(np.median([r[key] for r in rows]))

    print("\n=== aggregate over %d molecules (MMFF94 neutral ruler; strain = E above each ensemble's own min) ===" % len(rows))
    print(f"  native median energy - her median energy   : {med('d_med'):+.2f} kcal/mol  (>0 = native more strained)")
    print(f"  native p90 energy    - her p90 energy       : {med('d_p90'):+.2f} kcal/mol")
    print(f"  median STRAIN (E - ensemble min):  her {med('h_strmed'):.2f}   native {med('n_strmed'):.2f} kcal/mol")
    print(f"  p90    STRAIN (E - ensemble min):  her {med('h_strp90'):.2f}   native {med('n_strp90'):.2f} kcal/mol")
    print(f"  native ensemble min - her ensemble min      : {med('n_minus_h_min'):+.2f} kcal/mol  (<0 = native finds lower minimum)")
    fa = float(np.median([r["frac_above"] for r in rows]))
    fa_mean = float(np.mean([r["frac_above"] for r in rows]))
    print(f"  JUNK indicator: native confs above HER max energy: median {fa*100:.1f}%  mean {fa_mean*100:.1f}%")
    worst = sorted(rows, key=lambda r: -r["frac_above"])[:6]
    print("  most-junk molecules:", ", ".join(f"{r['mid']}({r['frac_above']*100:.0f}%)" for r in worst))
    nuff = sum(1 for r in rows if r["ff"] == "UFF")
    print(f"  (energy fn: MMFF for {len(rows)-nuff}, UFF fallback for {nuff})")
    json.dump(rows, open(HERE / "strain_compare.json", "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
