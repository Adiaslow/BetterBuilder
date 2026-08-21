//! Native Kekulization — resolve aromatic bonds to an alternating single/double (Kekulé) assignment,
//! RDKit-free. A faithful port of the chemistry in RDKit's `Kekulize.cpp`: `markDbondCands` decides
//! which aromatic atoms must take a double bond (the valence/charge rules that distinguish, e.g.,
//! pyridine-N, which does, from pyrrole-N, which does not), then a backtracking perfect matching
//! assigns one double bond to each such atom over the aromatic bonds.
//!
//! We deliberately do NOT reproduce RDKit's *canonical* Kekulé form (which would need canonical atom
//! ranking): the only consumer is SYBYL `aro6`, and that has been validated to be identical across
//! every Kekulé form (all 100 corpus molecules), so any valid form is correct. RDKit is used only to
//! validate this, never in the result.

use crate::smarts_match::Perceived;
use crate::valence::valence_list;

/// The primitives kekulization reads — a view usable before a full [`Perceived`] exists (RDKit runs
/// Kekulize mid-sanitize, before conjugation/hybridization). `total_valence` is explicit+implicit.
pub struct KekMol<'a> {
    pub atomic_numbers: &'a [u8],
    pub bonds: &'a [(usize, usize)],
    pub bond_order: &'a [u8],
    pub bond_aromatic: &'a [bool],
    pub aromatic: &'a [bool],
    pub charges: &'a [i8],
    pub total_valence: &'a [i32],
    /// Implicit + bracket hydrogens per atom. RDKit runs Kekulize pre-AddHs, so `markDbondCands` counts
    /// these toward an atom's bond order and degree; zero when called on a post-AddHs molecule.
    pub total_num_hs: &'a [i32],
}

/// Atomic numbers for which RDKit's `isEarlyAtom` is true (verbatim from `Atom.cpp`'s table). Used to
/// flip the formal-charge sign in the valence calculation.
const EARLY_ATOMS: &[u8] = &[
    3, 4, 5, 11, 12, 13, 19, 20, 21, 22, 30, 31, 32, 37, 38, 39, 40, 41, 48, 49, 50, 51, 55, 56, 57,
    58, 59, 60, 61, 72, 73, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 93, 104, 105, 106, 107, 108,
    109, 110, 111, 112, 113, 114, 115, 116, 117, 118,
];

fn is_early(z: u8) -> bool {
    EARLY_ATOMS.contains(&z)
}

/// Per-bond order after Kekulization: non-aromatic bonds keep their order; aromatic bonds become 1 or
/// 2. All atoms are assumed post-AddHs (no implicit hydrogens, no radicals) — the state this project
/// perceives in.
pub fn kekule_bond_orders(p: &Perceived) -> Vec<u8> {
    // a Perceived is post-AddHs: hydrogens are explicit neighbours already, so none to add here
    let no_implicit_hs = vec![0i32; p.atomic_numbers.len()];
    kekule_bond_orders_raw(&KekMol {
        atomic_numbers: &p.atomic_numbers,
        bonds: &p.bonds,
        bond_order: &p.bond_order,
        bond_aromatic: &p.bond_aromatic,
        aromatic: &p.aromatic,
        charges: &p.charges,
        total_valence: &p.total_valence,
        total_num_hs: &no_implicit_hs,
    })
}

/// Kekulize from primitives — see [`KekMol`].
pub fn kekule_bond_orders_raw(m: &KekMol) -> Vec<u8> {
    let n = m.atomic_numbers.len();
    let mut orders = m.bond_order.to_vec();
    // adjacency: per atom, (bond index, other atom)
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, &(a, b)) in m.bonds.iter().enumerate() {
        adj[a].push((bi, b));
        adj[b].push((bi, a));
    }
    // aromatic bonds start single
    for (o, &arom) in orders.iter_mut().zip(m.bond_aromatic) {
        if arom {
            *o = 1;
        }
    }
    let cand = mark_dbond_cands(m, &adj);
    let mut matched = vec![false; n];
    assign_matching(m, &adj, &cand, &mut matched, &mut orders);
    orders
}

/// `markDbondCands`: which aromatic atoms must take a double bond.
fn mark_dbond_cands(m: &KekMol, adj: &[Vec<(usize, usize)>]) -> Vec<bool> {
    let n = m.atomic_numbers.len();
    let mut cand = vec![false; n];
    for i in 0..n {
        if !m.aromatic[i] {
            continue;
        }
        let z = m.atomic_numbers[i];
        // sbo: aromatic bonds count 1 each; non-aromatic add their order; + totalNumHs
        let mut sbo = m.total_num_hs[i];
        let mut n_to_ignore = 0i32;
        for &(bi, _) in &adj[i] {
            if m.bond_aromatic[bi] {
                sbo += 1;
            } else {
                let c = i32::from(m.bond_order[bi]);
                sbo += c;
                if c == 0 {
                    n_to_ignore += 1;
                }
            }
        }
        let vals = valence_list(z).unwrap_or(&[]);
        let mut dv = i32::from(vals.first().copied().unwrap_or(0)); // default valence
        let mut chrg = i32::from(m.charges[i]);
        if is_early(z) {
            chrg = -chrg;
        }
        if z == 6 && chrg > 0 {
            chrg = -chrg;
        }
        dv += chrg;
        let tbo = m.total_valence[i];
        // bump dv through the valence list while total bond order exceeds it
        let mut vi = 1;
        while tbo > dv && vi < vals.len() && vals[vi] > 0 {
            dv = i32::from(vals[vi]) + chrg;
            vi += 1;
        }
        let total_degree = adj[i].len() as i32 - n_to_ignore + m.total_num_hs[i];
        // post-AddHs: no radicals; noImplicit is true (implicit H were made explicit)
        if total_degree >= dv {
            cand[i] = false;
        } else if dv == sbo + 1 {
            cand[i] = true;
        } else if dv == sbo + 2 {
            cand[i] = true; // noImplicit radical-allowing case
        }
    }
    cand
}

/// Assign one double bond to each candidate via a backtracking perfect matching over aromatic bonds.
fn assign_matching(
    m: &KekMol,
    adj: &[Vec<(usize, usize)>],
    cand: &[bool],
    matched: &mut [bool],
    orders: &mut [u8],
) -> bool {
    let next = (0..cand.len()).find(|&i| cand[i] && !matched[i]);
    let Some(i) = next else {
        return true; // every candidate matched
    };
    for &(bi, j) in &adj[i] {
        if m.bond_aromatic[bi] && cand[j] && !matched[j] {
            matched[i] = true;
            matched[j] = true;
            orders[bi] = 2;
            if assign_matching(m, adj, cand, matched, orders) {
                return true;
            }
            matched[i] = false;
            matched[j] = false;
            orders[bi] = 1;
        }
    }
    false
}
