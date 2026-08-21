"""Dump RDKit's kekulized bond orders and perceived aromaticity over a corpus.

    python3 dump_kekulize.py <corpus.smi> <out.json>

Two separable things, captured together because they come from one sanitization:

* `kekule_bond_orders` — the alternating single/double assignment `SANITIZE_KEKULIZE` produces.
  Aromatic bonds carry no integer order until this runs, and the assignment feeds both the bounds
  matrix (bond lengths) and the electron counting that aromaticity perception then does.
* `aromatic_atoms` / `aromatic_bonds` — what `SANITIZE_SETAROMATICITY` concludes, which is the
  perceived flags rather than what the SMILES happened to write.

Recorded per atom and per bond in post-`AddHs` index order, as sequences.
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ORDER = {
    Chem.BondType.SINGLE: 1,
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
    Chem.BondType.QUADRUPLE: 4,
    Chem.BondType.AROMATIC: 0,  # no integer order until kekulized
}


def perceive(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    molh = Chem.AddHs(mol)

    # perceived aromaticity, as sanitization left it
    aromatic_atoms = [int(a.GetIsAromatic()) for a in molh.GetAtoms()]
    aromatic_bonds = [int(b.GetIsAromatic()) for b in molh.GetBonds()]
    sanitized_orders = [ORDER.get(b.GetBondType(), -1) for b in molh.GetBonds()]

    # the kekulized assignment, on a copy so the flags above are untouched
    kek = Chem.Mol(molh)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    kekule_orders = [ORDER.get(b.GetBondType(), -1) for b in kek.GetBonds()]

    return {
        "n_atoms": molh.GetNumAtoms(),
        "n_bonds": molh.GetNumBonds(),
        "aromatic_atoms": aromatic_atoms,
        "aromatic_bonds": aromatic_bonds,
        "sanitized_bond_orders": sanitized_orders,
        "kekule_bond_orders": kekule_orders,
    }


def main(corpus, outpath):
    out, bad = {}, 0
    for line in open(corpus):
        f = line.split()
        if len(f) < 2:
            continue
        try:
            r = perceive(f[0])
        except Exception:
            r = None
        if r is None:
            bad += 1
            continue
        out[f[1]] = r
    with open(outpath, "w") as fh:
        json.dump(out, fh, sort_keys=True)
    arom_mols = sum(1 for v in out.values() if any(v["aromatic_atoms"]))
    arom_atoms = sum(sum(v["aromatic_atoms"]) for v in out.values())
    print(f"  {len(out)} molecules ({bad} failed), {arom_mols} with aromatic atoms, "
          f"{arom_atoms} aromatic atoms total")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
