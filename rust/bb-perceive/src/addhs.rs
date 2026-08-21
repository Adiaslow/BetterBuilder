//! Making hydrogens explicit, in RDKit's atom and bond order.
//!
//! This is what fixes the indexing every later stage means. RDKit appends hydrogens after all heavy
//! atoms, grouped by parent in parent order, bracket hydrogens before implicit ones, each with a
//! single bond added in that same sequence (`AddHs.cpp`). Both orders are reproduced here, because
//! both are read downstream: atom order by every per-atom array, and bond order by ring perception,
//! which walks an atom's neighbours in it.

use crate::hydrogens::Counts;
use crate::SmilesGraph;

/// A molecule with explicit hydrogens.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Explicit {
    /// Atomic numbers: the heavy atoms unchanged, then the appended hydrogens.
    pub atomic_numbers: Vec<u8>,
    /// Bonds in RDKit's order: the heavy-atom bonds, then one per appended hydrogen.
    pub bonds: Vec<(usize, usize)>,
    /// The heavy atom each appended hydrogen hangs off, in append order.
    pub h_parents: Vec<usize>,
    pub n_heavy: usize,
}

/// Append explicit hydrogens to a parsed graph.
pub fn add_hs(g: &SmilesGraph, counts: &Counts) -> Explicit {
    let n_heavy = g.atoms.len();
    let mut atomic_numbers: Vec<u8> = g.atoms.iter().map(|a| a.atomic_number).collect();
    let mut bonds = g.bonds.clone();
    let mut h_parents = Vec::new();

    for aidx in 0..g.atoms.len() {
        // append the total hydrogen count — bracket + implicit, after RDKit's `adjustHs` has moved any
        // implicit H lost to aromatization into the explicit count (all Hs on one parent are identical,
        // so the bracket-before-implicit ordering within a parent carries no downstream meaning).
        for _ in 0..counts.total_num_hs[aidx].max(0) {
            let new_idx = atomic_numbers.len();
            atomic_numbers.push(1);
            bonds.push((aidx, new_idx));
            h_parents.push(aidx);
        }
    }

    Explicit {
        atomic_numbers,
        bonds,
        h_parents,
        n_heavy,
    }
}
