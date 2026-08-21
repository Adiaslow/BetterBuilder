"""Dump RDKit's ring set as ordered walks, in ring order, over a corpus.

    python3 dump_ring_order.py <corpus.smi> <out.json>

`dump_rings.py` records the same rings as sorted sets, which is the right shape for the quantities
that only read membership. It is the wrong shape for `getExperimentalTorsions`, which walks
`RingInfo::atomRings()` positionally: the flat-ring torsions it emits are consecutive quadruples
`(i, i+1, i+2, i+3)` around each ring, and the first ring to claim a central bond locks it. Both the
order rings appear in and the direction each is walked therefore decide which quadruple is assigned.

Rings are recorded exactly as RDKit returns them — no sorting at either level. Atom indices are
post-`AddHs`, matching the indices torsions are reported in.
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def rings_for(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    molh = Chem.AddHs(mol)
    Chem.GetSymmSSSR(molh)
    return {
        "n_atoms": molh.GetNumAtoms(),
        # verbatim: ring order and walk direction are both load-bearing
        "rings": [list(r) for r in molh.GetRingInfo().AtomRings()],
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
    n_rings = sum(len(v["rings"]) for v in out.values())
    print(f"  {len(out)} molecules ({bad} unparseable), {n_rings} rings as ordered walks")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
