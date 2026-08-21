//! RDKit's `MolOps::cleanUp`: the charge-separation normalisations `sanitizeMol` runs first, before
//! the valence check or any perception. Three rewrites bring hypervalent neutral forms to the
//! zwitterionic forms RDKit's valence check accepts:
//!
//! - `nitrogensCleanup`: neutral 5-valent N with a double bond to O → `[N+]…[O-]` (nitro); and the
//!   azide tail `N=N#N` → `N=[N+]=[N-]`.
//! - `phosphorusCleanup`: neutral 5-valent, 3-coordinate P with one `=O` and one `=C`/`=N` →
//!   `[P+]…[O-]`.
//! - `halogenCleanup`: neutral Cl/Br/I of valence 3/5/7 whose neighbours are all O →
//!   `[X+n]([O-])…` (e.g. perchlorate `OCl(=O)(=O)=O` → `[Cl+3]([O-])([O-])([O-])O`).
//!
//! Operates on the as-written graph (aromaticity not yet perceived), mirroring RDKit's own order.

use crate::hydrogens::{contribution, explicit_valence};
use crate::{valence, SmilesGraph};
use yowl::feature::BondKind;

/// Non-metal elements — RDKit's `M` (metal) atom query is the complement of this set (Marvin's
/// definition, plus the dummy `#0`). Everything else counts as a metal.
const NON_METALS: [u8; 23] = [
    0, 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86,
];

/// RDKit's `QueryOps::isMetal`.
fn is_metal(z: u8) -> bool {
    !NON_METALS.contains(&z)
}

/// `noDative`: elements never made a dative donor (H, He, F, Ne).
fn no_dative(z: u8) -> bool {
    matches!(z, 1 | 2 | 9 | 10)
}

/// This atom's incident bonds as `(bond index, other atom index)`.
fn neighbors(g: &SmilesGraph, i: usize) -> Vec<(usize, usize)> {
    g.bonds
        .iter()
        .enumerate()
        .filter_map(|(bi, &(a, b))| {
            if a == i {
                Some((bi, b))
            } else if b == i {
                Some((bi, a))
            } else {
                None
            }
        })
        .collect()
}

/// Heavy-atom degree, matching RDKit's `getDegree()` on a graph whose hydrogens are implicit counts.
fn degree(g: &SmilesGraph, i: usize) -> usize {
    g.bonds.iter().filter(|&&(a, b)| a == i || b == i).count()
}

/// Whether bond `bi` is a single bond in RDKit's sense — order 1 and non-aromatic. An unwritten
/// (`Elided`) bond between two non-aromatic atoms is a single bond, which is how ligand→metal bonds
/// are usually written (e.g. ferrocene's ring-closure bonds to Fe).
fn is_single_bond(g: &SmilesGraph, bi: usize) -> bool {
    let (a, b) = g.bonds[bi];
    let both = g.atoms[a].aromatic_as_written && g.atoms[b].aromatic_as_written;
    (contribution(g.bond_kinds[bi], both) - 1.0).abs() < f64::EPSILON
}

/// Explicit valence of atom `i` from the graph's current bond/charge state — RDKit's
/// `calcExplicitValence(false)`, called at cleanUp time before aromaticity perception.
fn explicit_valence_of(g: &SmilesGraph, i: usize) -> i32 {
    let orders: Vec<f64> = neighbors(g, i)
        .iter()
        .map(|&(bi, o)| {
            if g.dative_donor[bi] == Some(i) {
                return 0.0; // this atom donates the dative bond: it contributes nothing here
            }
            let both = g.atoms[i].aromatic_as_written && g.atoms[o].aromatic_as_written;
            contribution(g.bond_kinds[bi], both)
        })
        .collect();
    let atom = &g.atoms[i];
    explicit_valence(
        atom.atomic_number,
        atom.charge,
        atom.bracket_hydrogens.unwrap_or(0),
        atom.aromatic_as_written,
        &orders,
    )
}

/// `nitrogensCleanup`: nitro and azide charge separation.
fn nitrogens_cleanup(g: &mut SmilesGraph) {
    // the neutral 5-valent N atoms RDKit considers, gathered first so the azide pass can revisit them
    let considered: Vec<usize> = (0..g.atoms.len())
        .filter(|&i| {
            g.atoms[i].atomic_number == 7
                && g.atoms[i].charge == 0
                && explicit_valence_of(g, i) == 5
        })
        .collect();
    for &i in &considered {
        // first neutral double-bonded O becomes [O-], the N becomes [N+]
        for (bi, o) in neighbors(g, i) {
            if g.atoms[o].atomic_number == 8
                && g.atoms[o].charge == 0
                && g.bond_kinds[bi] == BondKind::Double
            {
                g.bond_kinds[bi] = BondKind::Single;
                g.atoms[i].charge = 1;
                g.atoms[o].charge = -1;
                break;
            }
        }
    }
    // the azide tail: a triple bond to a neutral N drops to a double, N#N → [N+]=[N-]
    for &i in &considered {
        for (bi, o) in neighbors(g, i) {
            if g.atoms[o].atomic_number == 7
                && g.atoms[o].charge == 0
                && g.bond_kinds[bi] == BondKind::Triple
            {
                g.bond_kinds[bi] = BondKind::Double;
                g.atoms[i].charge = 1;
                g.atoms[o].charge = -1;
                break;
            }
        }
    }
}

/// `phosphorusCleanup`: neutral `C=P(=O)X` → `C=[P+]([O-])X`.
fn phosphorus_cleanup(g: &mut SmilesGraph, i: usize) {
    if g.atoms[i].charge != 0 || explicit_valence_of(g, i) != 5 || degree(g, i) != 3 {
        return;
    }
    let mut dbl_to_o: Option<(usize, usize)> = None;
    let mut has_double_to_c_or_n = false;
    for &(bi, o) in &neighbors(g, i) {
        let oz = g.atoms[o].atomic_number;
        if oz == 8 && g.atoms[o].charge == 0 && g.bond_kinds[bi] == BondKind::Double {
            dbl_to_o = Some((bi, o));
        } else if (oz == 6 || oz == 7) && degree(g, o) >= 2 && g.bond_kinds[bi] == BondKind::Double {
            has_double_to_c_or_n = true;
        }
    }
    if has_double_to_c_or_n {
        if let Some((bi, o)) = dbl_to_o {
            g.atoms[o].charge = -1;
            g.bond_kinds[bi] = BondKind::Single;
            g.atoms[i].charge = 1;
        }
    }
}

/// `halogenCleanup`: neutral all-oxygen Cl/Br/I of valence 3/5/7 → charge-separated oxyacid form.
fn halogen_cleanup(g: &mut SmilesGraph, i: usize) {
    let ev = explicit_valence_of(g, i);
    if g.atoms[i].charge != 0 || !(ev == 7 || ev == 5 || ev == 3) {
        return;
    }
    let nb = neighbors(g, i);
    if !nb.iter().all(|&(_, o)| g.atoms[o].atomic_number == 8) {
        return;
    }
    let mut formal_charge: i8 = 0;
    for (bi, o) in nb {
        if g.bond_kinds[bi] == BondKind::Double {
            g.bond_kinds[bi] = BondKind::Single;
            formal_charge += 1;
            g.atoms[o].charge = -1;
        }
    }
    g.atoms[i].charge = formal_charge;
}

/// `isHypervalentNonMetal`: a non-metal whose explicit valence exceeds the maximum allowed for its
/// charge-adjusted element — or, for cyclopentadienyl-type systems, equals it while aromatic and
/// 4-coordinate. `effAtomicNum = Z − charge` makes `N+` behave like `C`, etc.
fn is_hypervalent_non_metal(g: &SmilesGraph, i: usize) -> bool {
    let z = g.atoms[i].atomic_number;
    if is_metal(z) {
        return false;
    }
    let eff = i32::from(z) - i32::from(g.atoms[i].charge);
    if eff <= 0 {
        return false;
    }
    let max_v = match valence::valence_list(eff as u8) {
        Some(v) => i32::from(*v.last().unwrap_or(&-1)),
        None => -1,
    };
    if max_v <= 0 {
        return false;
    }
    let ev = explicit_valence_of(g, i);
    // getTotalDegree = heavy degree + bracket (explicit) H; the special case is for aromatic Cp-type
    // carbons, which carry their hydrogens explicitly when written aromatic.
    let total_degree = degree(g, i) + usize::from(g.atoms[i].bracket_hydrogens.unwrap_or(0));
    ev > max_v || (ev == max_v && g.atoms[i].aromatic_as_written && total_degree == 4)
}

/// Dative bonds incident on a metal atom (it is always the acceptor).
fn num_dative_bonds(g: &SmilesGraph, m: usize) -> usize {
    neighbors(g, m)
        .iter()
        .filter(|&&(bi, _)| g.dative_donor[bi].is_some())
        .count()
}

/// `metalBondCleanup`: retype one ligand→metal single bond from this hypervalent non-metal as
/// dative. Among candidate metals, RDKit picks the one with the fewest dative bonds (canonical rank
/// as the tie-breaker); we use the atom index as the rank surrogate (see [`clean_up_organometallics`]).
fn metal_bond_cleanup(g: &mut SmilesGraph, i: usize) {
    if !is_hypervalent_non_metal(g, i) || no_dative(g.atoms[i].atomic_number) {
        return;
    }
    let mut metals: Vec<(usize, usize)> = neighbors(g, i)
        .into_iter()
        .filter(|&(bi, o)| {
            is_single_bond(g, bi)
                && g.dative_donor[bi].is_none()
                && is_metal(g.atoms[o].atomic_number)
        })
        .collect();
    if metals.is_empty() {
        return;
    }
    // fewer dative bonds first; higher rank (index surrogate) breaks ties
    metals.sort_by(|&(_, m1), &(_, m2)| {
        num_dative_bonds(g, m1)
            .cmp(&num_dative_bonds(g, m2))
            .then(m2.cmp(&m1))
    });
    let (bi, _metal) = metals[0];
    g.dative_donor[bi] = Some(i);
}

/// `cleanUpOrganometallics`: retype ligand→metal single bonds as dative wherever a hypervalent
/// non-metal coordinates a metal, so the ligand atom's valence stops counting the metal bond (e.g.
/// ferrocene's cyclopentadienyl carbons). RDKit processes atoms in canonical-atom-rank order and
/// uses that rank to break ties between candidate metals; we substitute the atom index. The order
/// only changes the outcome when one non-metal is bonded to more than one metal by a single bond —
/// no molecule in the validation corpora does — so for all tested input this reproduces RDKit
/// exactly. A general bridged-polymetallic port would need RDKit's `Canon::rankMolAtoms`.
pub fn clean_up_organometallics(g: &mut SmilesGraph) {
    let needs_fixing = (0..g.atoms.len()).any(|i| {
        is_hypervalent_non_metal(g, i)
            && !no_dative(g.atoms[i].atomic_number)
            && neighbors(g, i)
                .iter()
                .any(|&(bi, o)| is_single_bond(g, bi) && is_metal(g.atoms[o].atomic_number))
    });
    if !needs_fixing {
        return;
    }
    for i in 0..g.atoms.len() {
        metal_bond_cleanup(g, i);
    }
}

/// Run RDKit's `cleanUp` in its order: nitrogens first, then per-atom phosphorus/halogen rewrites.
pub fn clean_up(g: &mut SmilesGraph) {
    nitrogens_cleanup(g);
    for i in 0..g.atoms.len() {
        match g.atoms[i].atomic_number {
            15 => phosphorus_cleanup(g, i),
            17 | 35 | 53 => halogen_cleanup(g, i),
            _ => {}
        }
    }
}
