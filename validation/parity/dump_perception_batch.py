"""Dump RDKit's perception for many molecules, for comparison against bb-rdkit's MoleculeSpec.

Runs inside the pipeline container, so the RDKit is the same patched build the oracle uses.
Every value comes from RDKit; nothing is reimplemented.

    python3 dump_perception_batch.py <corpus.smi> <out.json> [count]
"""
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import check  # single source for mol_id (the corpus key)
from rdkit import Chem
from rdkit.Chem import rdDistGeom


def perceive(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    n = mol.GetNumAtoms()
    bounds = rdDistGeom.GetMoleculeBoundsMatrix(mol, useMacrocycle14config=True)
    ps = rdDistGeom.ETKDGv3()
    ps.useMacrocycleTorsions = True
    ps.useExpTorsionAnglePrefs = True
    tors = rdDistGeom.GetExperimentalTorsions(mol, ps) or []
    return {
        "n_atoms": n,
        "formal_charge": Chem.rdmolops.GetFormalCharge(mol),
        "atomic_numbers": [a.GetAtomicNum() for a in mol.GetAtoms()],
        "bonds": sorted(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()),
        "bounds": [[float(bounds[i][j]) for j in range(n)] for i in range(n)],
        "torsions": sorted(
            [sorted(t["atomIndices"]) + [round(float(v), 6) for v in t["V"]] for t in tors]
        ),
    }


def main(smipath, outpath, count=None):
    out = {}
    for line in open(smipath):
        f = line.split()
        if len(f) < 2:
            continue
        out[check.mol_id(f[1])] = perceive(f[0])
        if count and len(out) >= int(count):
            break
    with open(outpath, "w") as fh:
        json.dump(out, fh)
    print(f"  perceived {len(out)} molecules")


if __name__ == "__main__":
    main(*sys.argv[1:4])
