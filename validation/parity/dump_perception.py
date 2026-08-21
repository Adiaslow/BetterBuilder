"""Dump RDKit's perception of a molecule, for comparison against bb-rdkit's MoleculeSpec.

Runs inside the pipeline container, so the RDKit used here is the same patched build the oracle
uses. Every value comes from RDKit itself; nothing is reimplemented.

    python3 dump_perception.py "<smiles>" <out.json>

bb-rdkit extracts these same quantities through the cxx bridge. Comparing them checks the
extraction — the foundation every later stage rests on.
"""
import json
import sys

from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdMolDescriptors


def main(smiles, outpath):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))

    # Bounds matrix, exactly as the embedder builds it (ETKDGv3 uses macrocycle 1-4 bounds).
    bounds = rdDistGeom.GetMoleculeBoundsMatrix(mol, useMacrocycle14config=True)
    n = mol.GetNumAtoms()

    # CrystalFF experimental torsions, carrying the amide patch.
    ps = rdDistGeom.ETKDGv3()
    ps.useMacrocycleTorsions = True
    ps.useExpTorsionAnglePrefs = True
    tors = rdDistGeom.GetExperimentalTorsions(mol, ps)

    out = {
        "smiles": smiles,
        "n_atoms": n,
        "formal_charge": Chem.rdmolops.GetFormalCharge(mol),
        "atomic_numbers": [a.GetAtomicNum() for a in mol.GetAtoms()],
        "bonds": sorted(
            [sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()]
        ),
        # row-major n*n; RDKit's convention is upper = upper bound, lower = lower bound
        "bounds": [[float(bounds[i][j]) for j in range(n)] for i in range(n)],
        "experimental_torsions": [
            {
                "smarts": t["smarts"],
                "atoms": list(t.get("atomIndices", t.get("atoms", []))),
                "V": [float(x) for x in t["V"]],
                "signs": [int(x) for x in t["signs"]],
            }
            for t in (tors if tors is not None else [])
        ],
        "n_rotatable": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "ring_info": [sorted(r) for r in Chem.GetSymmSSSR(mol)],
    }
    with open(outpath, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(f"  {n} atoms, charge {out['formal_charge']}, "
          f"{len(out['experimental_torsions'])} torsions, {len(out['ring_info'])} rings")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
