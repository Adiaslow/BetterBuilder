"""Dump RDKit's molecule after `AddHs`, over a corpus.

    python3 dump_addhs.py <corpus.smi> <out.json>

`AddHs` fixes the atom ordering every downstream index means: the bounds matrix is n x n over it,
torsions and chiral sets are atom quadruples in it, and AMSOL reports its per-atom charges back in
it. An off-by-one here is silent and misaligns everything after.

Both atoms and bonds are recorded **in RDKit's order**, not sorted. Bond order is load-bearing —
ring perception walks an atom's neighbours in bond order, and among equal-size rings that decides
which is found — so a list that matches as a set but not as a sequence is not a match.
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def perceive(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    n_heavy = mol.GetNumAtoms()
    molh = Chem.AddHs(mol)
    return {
        "n_heavy": n_heavy,
        "n_atoms": molh.GetNumAtoms(),
        "z": [a.GetAtomicNum() for a in molh.GetAtoms()],
        # sequence, not a set: index i is RDKit's bond i
        "bonds": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in molh.GetBonds()],
        # the heavy atom each appended hydrogen hangs off, in append order
        "h_parents": [
            n.GetIdx()
            for a in molh.GetAtoms()
            if a.GetAtomicNum() == 1
            for n in a.GetNeighbors()
        ],
    }


def main(corpus, outpath):
    out, bad = {}, 0
    for line in open(corpus):
        f = line.split()
        if len(f) < 2:
            continue
        r = perceive(f[0])
        if r is None:
            bad += 1
            continue
        out[f[1]] = r
    with open(outpath, "w") as fh:
        json.dump(out, fh, sort_keys=True)
    atoms = sum(v["n_atoms"] for v in out.values())
    bonds = sum(len(v["bonds"]) for v in out.values())
    print(f"  {len(out)} molecules ({bad} unparseable), {atoms} atoms, {bonds} bonds after AddHs")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
