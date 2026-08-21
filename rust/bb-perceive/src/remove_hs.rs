//! RDKit's `MolOps::removeHs`, as `MolFromSmiles` runs it (default params + `updateExplicitCount`).
//!
//! A hydrogen written explicitly in a SMILES (`[C@](...)[H]`) is folded back into its neighbour's
//! explicit-H count and dropped from the graph; `bb_perceive::addhs` re-adds it later, so the H ends
//! up at RDKit's canonical (post-AddHs) index rather than its written one. When the neighbour is a
//! stereocentre, removing the H reorders the bonds about it, which can flip the chiral tag — RDKit
//! moves the H's bond to last and inverts the tag if that permutation is odd, and so do we.

use crate::SmilesGraph;
use yowl::feature::BondKind;

/// Whether hydrogen atom `h` is one `removeHs` drops, matching `shouldRemoveH` with the SMILES
/// parser's parameters (degree 1, not a hydride, neighbour not a dummy, bond not stereo-defining).
fn removable(g: &SmilesGraph, h: usize, incident: &[Vec<usize>]) -> bool {
    if g.atoms[h].atomic_number != 1 {
        return false;
    }
    if g.atoms[h].charge == -1 {
        return false; // removeHydrides = false: keep [H-]
    }
    if g.atoms[h].isotope {
        return false; // removeIsotopes = false: keep deuterium/tritium
    }
    if incident[h].len() != 1 {
        return false; // removeDegreeZero / removeHigherDegrees = false: only degree-1 H
    }
    let bi = incident[h][0];
    // removeDefiningBondStereo = false: keep an H whose bond carries a `/`,`\` direction
    if matches!(g.bond_kinds[bi], BondKind::Up | BondKind::Down) {
        return false;
    }
    let (a, b) = g.bonds[bi];
    let other = if a == h { b } else { a };
    // neighbour must be a real heavy atom: Z > 1 excludes both a dummy (removeDummyNeighbors=false)
    // and another hydrogen (removeOnlyHNeighbors=false).
    g.atoms[other].atomic_number > 1
}

/// Fold explicit hydrogens back into their neighbours (RDKit's `removeHs`), returning the reduced
/// graph. Atom indices are compacted; `chirality` and `dative_donor` are remapped to match.
///
/// `n_chain` is the number of leading chain bonds in `g.bonds` (ring-closure bonds are appended
/// after them). The chirality perturbation counts only chain bonds written after the H — a ring
/// closure is written at the atom's token, i.e. before any branch or explicit H, so it never sits
/// "after" the H even though `g.bonds` stores it later.
pub fn remove_hs(g: SmilesGraph, n_chain: usize) -> SmilesGraph {
    let n = g.atoms.len();
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        incident[a].push(bi);
        incident[b].push(bi);
    }

    let remove: Vec<bool> = (0..n).map(|h| removable(&g, h, &incident)).collect();
    if !remove.iter().any(|&r| r) {
        return g;
    }

    let mut chirality = g.chirality.clone();
    let mut atoms = g.atoms.clone();
    // 1. account for each removed H at its neighbour: bump the explicit-H count, and — if the
    //    neighbour is a stereocentre — invert its tag when moving the H's bond to last is an odd
    //    permutation of the neighbour's bonds (RDKit's `getPerturbationOrder`).
    for h in 0..n {
        if !remove[h] {
            continue;
        }
        let bi = incident[h][0];
        let (a, b) = g.bonds[bi];
        let p = if a == h { b } else { a };
        atoms[p].bracket_hydrogens = Some(atoms[p].bracket_hydrogens.unwrap_or(0) + 1);
        chirality.bump_num_hs(p);
        if chirality.is_chiral(p) {
            // swaps to move the H bond to last, in parse order: count P's chain bonds written after
            // the H (ring-closure bonds are written before, so they don't count).
            let swaps = incident[p]
                .iter()
                .filter(|&&x| x != bi && x < n_chain && x > bi)
                .count();
            if swaps % 2 == 1 {
                chirality.invert(p);
            }
        }
    }

    // 2. compact the surviving atoms and build the old→new index map
    let mut new_index = vec![usize::MAX; n];
    let mut keep: Vec<usize> = Vec::new();
    for i in 0..n {
        if !remove[i] {
            new_index[i] = keep.len();
            keep.push(i);
        }
    }
    let atoms: Vec<_> = keep.iter().map(|&i| atoms[i].clone()).collect();

    // 3. rebuild bonds (dropping the removed H bonds) and remap endpoints; carry kinds + dative
    let mut bonds = Vec::new();
    let mut bond_kinds = Vec::new();
    let mut dative_donor = Vec::new();
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        if remove[a] || remove[b] {
            continue;
        }
        bonds.push((new_index[a], new_index[b]));
        bond_kinds.push(g.bond_kinds[bi]);
        dative_donor.push(g.dative_donor[bi].map(|d| new_index[d]));
    }

    chirality.compact(&new_index, keep.len());

    SmilesGraph {
        atoms,
        bonds,
        bond_kinds,
        dative_donor,
        chirality,
    }
}
