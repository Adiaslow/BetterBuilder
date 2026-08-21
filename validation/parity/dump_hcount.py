"""Dump RDKit's per-atom valence and hydrogen counts over a corpus.

    python3 dump_hcount.py <corpus.smi> <out.json>

These are the values `SANITIZE_PROPERTIES` (`updatePropertyCache`) computes, which run before ring
perception and decide how many hydrogens `AddHs` appends — and therefore the atom ordering every
downstream index rides on.

Values are per atom in pre-AddHs order, so they line up with what `bb-perceive` produces from the
SMILES. Recorded as sequences, not totals: a molecule whose hydrogen count is right in aggregate but
wrong per atom produces a different AddHs ordering, and every later comparison would be misaligned.
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def perceive(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "n_atoms": mol.GetNumAtoms(),
        "z": [a.GetAtomicNum() for a in mol.GetAtoms()],
        "charge": [a.GetFormalCharge() for a in mol.GetAtoms()],
        "degree": [a.GetDegree() for a in mol.GetAtoms()],
        "aromatic": [int(a.GetIsAromatic()) for a in mol.GetAtoms()],
        "explicit_valence": [a.GetExplicitValence() for a in mol.GetAtoms()],
        "implicit_valence": [a.GetImplicitValence() for a in mol.GetAtoms()],
        "num_explicit_hs": [a.GetNumExplicitHs() for a in mol.GetAtoms()],
        "total_num_hs": [a.GetTotalNumHs() for a in mol.GetAtoms()],
        "radicals": [a.GetNumRadicalElectrons() for a in mol.GetAtoms()],
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
    hs = sum(sum(v["total_num_hs"]) for v in out.values())
    print(f"  {len(out)} molecules ({bad} unparseable), {atoms} heavy atoms, {hs} hydrogens")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
