"""Dump RDKit's conjugation and hybridization over a corpus.

    python3 dump_hybrid.py <corpus.smi> <out.json>

Hybridization is needed well beyond its obvious uses: UFF atom typing keys on it, so it reaches the
bond lengths and angles in `setTopolBounds`, and the UFF inversion terms only apply to SP2 centres.
It depends in turn on conjugation, which `setConjugation` marks per bond.

Per atom and per bond in post-AddHs order, as sequences.
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HYB = {
    Chem.HybridizationType.UNSPECIFIED: 0, Chem.HybridizationType.S: 1,
    Chem.HybridizationType.SP: 2, Chem.HybridizationType.SP2: 3,
    Chem.HybridizationType.SP3: 4, Chem.HybridizationType.SP3D: 5,
    Chem.HybridizationType.SP3D2: 6,
}


def perceive(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    molh = Chem.AddHs(mol)
    return {
        "n_atoms": molh.GetNumAtoms(),
        "hybridization": [HYB.get(a.GetHybridization(), -1) for a in molh.GetAtoms()],
        "conjugated": [int(b.GetIsConjugated()) for b in molh.GetBonds()],
        "total_degree": [a.GetTotalDegree() for a in molh.GetAtoms()],
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
    from collections import Counter
    c = Counter(h for v in out.values() for h in v["hybridization"])
    conj = sum(sum(v["conjugated"]) for v in out.values())
    print(f"  {len(out)} molecules ({bad} failed), {conj} conjugated bonds")
    print(f"  hybridization histogram (0=unspec,2=SP,3=SP2,4=SP3): {dict(sorted(c.items()))}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
