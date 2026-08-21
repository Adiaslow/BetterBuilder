#!/usr/bin/env python3
# Divergence-hunting generator (goal criterion 5: "inputs are divergence-hunted with measured
# coverage, so 'no diffs' means 'searched hard,' not 'the sample passed'"). CSmith-style structure-
# aware mutation: start from her REAL macrocycle workload (seeds_5000) and apply small, valence-valid,
# macrocycle-preserving mutations to probe the perception/FF/bounds code near — but off — the sampled
# corpus. The point is adversarial COVERAGE: exercise element/bond/stereo/ring combinations the fixed
# sample never contained, then diff native vs the RDKit bridge on every one (differential.py).
#
# Deterministic given SEED (no RNG-in-transcript concerns): a fixed-seed PRNG drives mutation choices.
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/hunt/mutate_corpus.py <seeds.smi> <out.smi> [n_per_seed] [rng_seed]
import collections
import random
import statistics
import sys
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

# element -> (atomic num, typical max valence) for structure-preserving swaps
SWAP = {6: 4, 7: 3, 8: 2, 16: 2}
ADD_ATOMS = [6, 7, 8, 9, 17]   # C N O F Cl as new substituents
MIN_MACRO = 9                  # keep at least one macrocyclic ring (>=9) so we stay in-domain


def max_ring(mol):
    ri = mol.GetRingInfo()
    return max((len(r) for r in ri.AtomRings()), default=0)


def valid_macrocycle(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return smi if max_ring(m) >= MIN_MACRO else None


def free_valence(atom):
    return atom.GetTotalNumHs()


def mutate(mol, rng):
    """Apply ONE random valence-valid, macrocycle-preserving mutation; return SMILES or None."""
    m = Chem.RWMol(mol)
    op = rng.choice(["swap", "add", "del", "bond"])
    atoms = list(m.GetAtoms())
    try:
        if op == "swap":
            a = rng.choice(atoms)
            cur = a.GetAtomicNum()
            cands = [z for z in SWAP if z != cur and a.GetDegree() <= SWAP[z]]
            if not cands:
                return None
            a.SetAtomicNum(rng.choice(cands))
            a.SetNumExplicitHs(0)
            a.SetNoImplicit(False)
        elif op == "add":
            host = rng.choice([a for a in atoms if free_valence(a) > 0] or atoms)
            if free_valence(host) <= 0:
                return None
            z = rng.choice(ADD_ATOMS)
            idx = m.AddAtom(Chem.Atom(z))
            m.AddBond(host.GetIdx(), idx, Chem.BondType.SINGLE)
        elif op == "del":
            terminals = [a for a in atoms if a.GetDegree() == 1 and not a.GetIsAromatic()]
            if not terminals:
                return None
            m.RemoveAtom(rng.choice(terminals).GetIdx())
        elif op == "bond":
            b = rng.choice(list(m.GetBonds()))
            if b.GetIsAromatic() or b.IsInRing():
                return None  # don't disturb ring/aromatic bonds (keep the macrocycle intact)
            i, j = b.GetBeginAtom(), b.GetEndAtom()
            if b.GetBondType() == Chem.BondType.SINGLE and free_valence(i) > 0 and free_valence(j) > 0:
                b.SetBondType(Chem.BondType.DOUBLE)
            elif b.GetBondType() == Chem.BondType.DOUBLE:
                b.SetBondType(Chem.BondType.SINGLE)
            else:
                return None
        m2 = m.GetMol()
        Chem.SanitizeMol(m2)
        return valid_macrocycle(Chem.MolToSmiles(m2)), op
    except (Chem.AtomValenceException, Chem.KekulizeException, ValueError, RuntimeError):
        return None


def main():
    seeds_path = sys.argv[1] if len(sys.argv) > 1 else "validation/seeds_5000.smi"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "validation/hunt/generated.smi"
    n_per = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 20210185)

    seeds = []
    for line in open(seeds_path):
        t = line.split()
        if t:
            seeds.append(t[0])
    seen = set(Chem.CanonSmiles(s) for s in seeds if Chem.MolFromSmiles(s))
    generated, op_counts = [], collections.Counter()
    tries = 0
    for smi in seeds:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        for _ in range(n_per):
            tries += 1
            r = mutate(mol, rng)
            if not r or r[0] is None:
                continue
            cand, op = r
            can = Chem.CanonSmiles(cand)
            if can in seen:
                continue
            seen.add(can)
            generated.append(can)
            op_counts[op] += 1

    with open(out_path, "w") as w:
        for i, s in enumerate(generated):
            w.write(f"{s} gen{i:06d}\n")

    # measured coverage: how hard did we search, and over what structural range?
    ringsz, hetero, natoms = collections.Counter(), collections.Counter(), []
    for s in generated:
        m = Chem.MolFromSmiles(s)
        ringsz[max_ring(m)] += 1
        hetero[sum(1 for a in m.GetAtoms() if a.GetAtomicNum() not in (1, 6))] += 1
        natoms.append(m.GetNumAtoms())
    print(f"generated {len(generated)} unique valid macrocycles from {len(seeds)} seeds "
          f"({tries} mutation attempts, {100.0*len(generated)/max(tries,1):.0f}% yield)")
    print(f"  mutation ops: {dict(op_counts)}")
    print(f"  macro ring-size coverage: {dict(sorted(ringsz.items()))}")
    print(f"  heteroatom-count coverage: {dict(sorted(hetero.items()))}")
    print(f"  atom-count range: {min(natoms)}..{max(natoms)} (median {statistics.median(natoms):.0f})")
    print(f"  -> {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
