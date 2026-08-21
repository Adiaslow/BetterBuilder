#!/usr/bin/env python
"""Dump the mol2db2 stage's input AND every hierarchy internal to JSON, using the container's own
authoritative code. This is the gate oracle for the Rust port: each stage (rigid structures, position
counts, rigid component, confs, sets, output coords, broken sets) is compared against these arrays,
not just the final db2.

    python dump_mol2db2.py <mol2> <solv> <out.json>

Run inside the pipeline container ($DOCKBASE/ligand/mol2db2_py3_strain on sys.path)."""
import json
import sys

import mol2 as mol2mod
import solv as solvmod
import clash as clashmod
import hierarchy as hiermod

DOCKBASE = "/soft/DOCK-macro-latest6"

def main(mol2file, solvfile, outfile):
    m = mol2mod.Mol2(mol2file, nameFileName=None)
    m.convertDockTypes(None)
    m.addColors(None)
    s = solvmod.Solv(solvfile)
    clashDecider = clashmod.Clash(DOCKBASE + "/ligand/mol2db2/clashfile.txt")
    h = hiermod.Hierarchy(m, clashDecider, tolerance=0.001, solvdata=s)

    out = {
        "input": {
            "name": m.name, "protName": m.protName, "smiles": m.smiles, "longname": m.longname,
            "atomNum": m.atomNum, "atomName": m.atomName, "atomType": m.atomType,
            "dockNum": m.dockNum, "colorNum": m.colorNum,
            "bondNum": m.bondNum, "bondStart": m.bondStart, "bondEnd": m.bondEnd, "bondType": m.bondType,
            "atomXyz": m.atomXyz,
            "inputTotalStrain": m.inputTotalStrain, "inputMaxStrain": m.inputMaxStrain,
            "inputHydrogens": m.inputHydrogens,
        },
        "solv": {
            "totalCharge": s.totalCharge, "totalPolarSolv": s.totalPolarSolv,
            "totalApolarSolv": s.totalApolarSolv, "totalSolv": s.totalSolv,
            "totalSurface": s.totalSurface,
            "charge": s.charge, "polarSolv": s.polarSolv, "apolarSolv": s.apolarSolv,
            "solv": s.solv, "surface": s.surface,
        },
        "hier": {
            "rigidStructures": h.rigidStructures, "rigidStructureCount": h.rigidStructureCount,
            "posCount": h.posCount, "outAtoms": h.outAtoms, "numConfs": h.numConfs,
            "heavyRigidCount": h.heavyRigidCount, "heavyRigidAtomNums": list(h.heavyRigidAtomNums),
            "outAtomOrigAtom": h.outAtomOrigAtom, "outAtomConfNum": h.outAtomConfNum,
            "outAtomXYZ": h.outAtomXYZ,
            "confNumAtomList": h.confNumAtomList,
            "setToConfs": {str(k): v for k, v in h.setToConfs.items()},
            "brokenSets": sorted(h.brokenSets),
        },
    }
    with open(outfile, "w") as f:
        json.dump(out, f)
    print(f"dumped {outfile}: {len(m.atomNum)} atoms, {len(m.atomXyz)} confs, "
          f"{h.numConfs} hierConfs, {h.rigidStructureCount} rigidStructs, {len(h.brokenSets)} broken")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
