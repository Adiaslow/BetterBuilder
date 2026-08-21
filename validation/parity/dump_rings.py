"""Dump RDKit's symmetrized SSSR for a corpus, for comparison against bb-perceive.

    python3 dump_rings.py <corpus.smi> <out.json>

`Chem.GetSymmSSSR` is what `SANITIZE_SYMMRINGS` stores on the molecule, so it is the ring set every
downstream query reads: `IsInRing`, ring-size SMARTS such as `r{9-}`, `count_exo_rotatable`'s ring
atoms, and the largest-ring pin selection.

Rings are emitted as sorted atom-index lists, themselves sorted, so the comparison is set equality
and does not depend on the order RDKit happens to return them in. Atom indices are pre-AddHs, which
is the order `bb-perceive` produces from the SMILES.
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def rings_for(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rings = [sorted(int(i) for i in r) for r in Chem.GetSymmSSSR(mol)]
    rings.sort()
    return {
        "n_atoms": mol.GetNumAtoms(),
        "n_bonds": mol.GetNumBonds(),
        "rings": rings,
    }


def main(corpus, outpath):
    out, bad = {}, 0
    for line in open(corpus):
        f = line.split()
        if len(f) < 2:
            continue
        r = rings_for(f[0])
        if r is None:
            bad += 1
            continue
        out[f[1]] = r
    with open(outpath, "w") as fh:
        json.dump(out, fh, sort_keys=True)
    sizes = {}
    for v in out.values():
        for ring in v["rings"]:
            sizes[len(ring)] = sizes.get(len(ring), 0) + 1
    print(f"  {len(out)} molecules ({bad} unparseable), "
          f"{sum(len(v['rings']) for v in out.values())} rings")
    print(f"  ring-size histogram: {dict(sorted(sizes.items()))}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
