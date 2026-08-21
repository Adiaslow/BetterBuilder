//! The two-stage recipe parameters and the pinned core.
//!
//! `count_exo_rotatable` is her function (`build_ligands.py`): single bonds that are not themselves
//! ring bonds but touch a ring atom. It runs on the molecule as parsed from SMILES, before
//! hydrogens are made explicit, and reads the symmetrized SSSR — so it needs ring perception to be
//! right, which is why it lands here rather than in the bridge.
//!
//! The pinned core is the largest ring: her seed loop fixes those atoms and samples the rest.

use crate::hydrogens::contribution;
use crate::SmilesGraph;

/// Single bonds outside a ring that touch a ring atom.
pub fn count_exo_rotatable(g: &SmilesGraph, rings: &[Vec<usize>]) -> usize {
    let ring_atoms: std::collections::HashSet<usize> =
        rings.iter().flat_map(|r| r.iter().copied()).collect();
    let in_ring = crate::sssr::ring_bond_flags_masked(g.atoms.len(), &g.bonds, &g.dative_donor);

    let mut exo = 0;
    for (bi, &(i, j)) in g.bonds.iter().enumerate() {
        if in_ring[bi] {
            continue;
        }
        let both_aromatic = g.atoms[i].aromatic_as_written && g.atoms[j].aromatic_as_written;
        // single bonds only; an aromatic bond contributes 1.5 and is excluded
        if (contribution(g.bond_kinds[bi], both_aromatic) - 1.0).abs() > f64::EPSILON {
            continue;
        }
        if ring_atoms.contains(&i) || ring_atoms.contains(&j) {
            exo += 1;
        }
    }
    exo
}

/// Core seeds and sidechain conformers, from her thresholds.
pub fn recipe(exo: usize) -> (u32, u32) {
    if exo <= 2 {
        (20, 10)
    } else {
        (10, 20)
    }
}

/// The pinned core: atoms of the largest ring.
pub fn core_atoms(rings: &[Vec<usize>]) -> Vec<usize> {
    // the FIRST largest ring, matching Python's max(sssr, key=len)
    let mut core = rings
        .iter()
        .fold(None::<&Vec<usize>>, |best, r| match best {
            Some(b) if b.len() >= r.len() => Some(b),
            _ => Some(r),
        })
        .cloned()
        .unwrap_or_default();
    core.sort_unstable();
    core
}

/// The full core-pin recipe set: the largest ring's atoms plus exocyclic pins — an exocyclic `=O`
/// (carbonyl) on any core atom, and an exocyclic C bonded to a core ring N. Reproduces
/// `build_ligands.py`'s cmap construction; returned sorted (RDKit iterates the pin set as a
/// `std::set`).
pub fn pin_atoms(g: &SmilesGraph, rings: &[Vec<usize>], bond_aromatic: &[bool]) -> Vec<usize> {
    use yowl::feature::BondKind;
    let core = core_atoms(rings);
    if core.is_empty() {
        return Vec::new();
    }
    let core_set: std::collections::HashSet<usize> = core.iter().copied().collect();
    let n = g.atoms.len();
    // adjacency with a genuine-double-bond flag per incident bond. RDKit's pin uses
    // `getBondTypeAsDouble() == 2.0`, which is 1.5 for an AROMATIC bond — so a PERCEIVED-aromatic bond
    // (e.g. a pyrylium C=O⁺ written in Kekulé form) is not a double here, and its O is not pinned.
    let mut adj: Vec<Vec<(usize, bool)>> = vec![Vec::new(); n];
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        let is_double = g.bond_kinds[bi] == BondKind::Double && !bond_aromatic[bi];
        adj[a].push((b, is_double));
        adj[b].push((a, is_double));
    }
    let mut pin: std::collections::BTreeSet<usize> = core.iter().copied().collect();
    for &k in &core {
        for &(nb, is_double) in &adj[k] {
            // exocyclic double-bonded O (e.g. a carbonyl)
            if g.atoms[nb].atomic_number == 8 && is_double {
                pin.insert(nb);
            }
            // exocyclic C bonded to a core ring N
            if g.atoms[k].atomic_number == 7
                && g.atoms[nb].atomic_number == 6
                && !core_set.contains(&nb)
            {
                pin.insert(nb);
            }
        }
    }
    pin.into_iter().collect()
}
