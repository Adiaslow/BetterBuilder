//! Explicit and implicit valence, and the hydrogen count that follows from them.
//!
//! These are the values `SANITIZE_PROPERTIES` computes before ring perception, and they decide how
//! many hydrogens `AddHs` appends — hence the atom ordering every downstream index rides on.
//!
//! Two details carry most of the behaviour. An atom's valence list is consulted for the element
//! shifted by its formal charge (`atomic number - charge`), which is why quaternary N+ takes four
//! bonds: it reads carbon's list. And phosphorus/sulfur (and As/Se) get an exemption from that
//! shift when it would push them past their own hypervalent range, so a sulfonamide S reaches 6.

use crate::valence;

/// Bond order as RDKit accumulates it; aromatic bonds contribute 1.5.
pub const AROMATIC_CONTRIB: f64 = 1.5;

/// RDKit clamps the charge-shifted element to the table's range.
fn effective_atomic_num(z: u8, charge: i8) -> u8 {
    let eff = i32::from(z) - i32::from(charge);
    eff.clamp(0, 118) as u8
}

/// Whether the charge shift is skipped so the element keeps its own hypervalent range.
fn can_be_hypervalent(z: u8, effective: u8) -> bool {
    (effective > 16 && (z == 15 || z == 16)) || (effective > 34 && (z == 33 || z == 34))
}

/// Does this element's valence list constrain it at all? A sole `-1` means unconstrained.
fn is_constrained(z: u8) -> bool {
    match valence::valence_list(z) {
        Some(v) => v.len() > 1 || v.first() != Some(&-1),
        None => false,
    }
}

fn default_valence(z: u8) -> i32 {
    valence::valence_list(z)
        .and_then(|v| v.first())
        .map_or(-1, |&v| i32::from(v))
}

/// RDKit's strict valence check (`calculateExplicitValence`, `strict=true`): reject an atom whose
/// rounded explicit valence exceeds the maximum its element allows. This is the sanitization step
/// that makes `MolFromSmiles` return null for `[CH5]`, `[ClH2]`, pentavalent carbon, etc. — perception
/// must refuse exactly what RDKit refuses. Returns the offending atom index, or `Ok` if all valid.
pub fn check_valences(g: &crate::SmilesGraph) -> Result<(), usize> {
    let mut orders: Vec<Vec<f64>> = vec![Vec::new(); g.atoms.len()];
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        let both = g.atoms[a].aromatic_as_written && g.atoms[b].aromatic_as_written;
        let c = contribution(g.bond_kinds[bi], both);
        // a dative bond contributes nothing to the donor (RDKit's `Bond::DATIVE` valence contrib)
        orders[a].push(if g.dative_donor[bi] == Some(a) { 0.0 } else { c });
        orders[b].push(if g.dative_donor[bi] == Some(b) { 0.0 } else { c });
    }
    for (i, atom) in g.atoms.iter().enumerate() {
        let z = atom.atomic_number;
        let chg = atom.charge;
        // an element with no valence list (`-1` sentinel, e.g. metals) is never rejected on valence
        let ovalens = match valence::valence_list(z) {
            Some(v) => v,
            None => continue,
        };
        let ovback = i32::from(*ovalens.last().unwrap_or(&-1));
        let bracket = atom.bracket_hydrogens.unwrap_or(0);
        let res = explicit_valence(z, chg, bracket, atom.aromatic_as_written, &orders[i]);
        let eff = if is_constrained(z) { effective_atomic_num(z, chg) } else { z };
        let valens = valence::valence_list(eff).unwrap_or(ovalens);
        let mut max_valence = i32::from(*valens.last().unwrap_or(&-1));
        let mut offset = 0i32;
        // negatively-charged P/S/As/Se keep their own (element) hypervalent range and shift by charge
        if can_be_hypervalent(z, eff) {
            max_valence = ovback;
            offset -= i32::from(chg);
        }
        // RDKit historically accepts two-coordinate [H-]
        if z == 1 && chg == -1 {
            max_valence = 2;
        }
        if max_valence >= 0 && ovback >= 0 && (res + offset) > max_valence {
            return Err(i);
        }
    }
    Ok(())
}

/// Explicit valence: bond orders plus bracket hydrogens, with RDKit's aromatic correction.
///
/// `bond_orders` are the valence contributions of this atom's bonds (1.5 for aromatic).
pub fn explicit_valence(
    z: u8,
    charge: i8,
    bracket_hydrogens: u8,
    aromatic: bool,
    bond_orders: &[f64],
) -> i32 {
    let mut accum: f64 = bond_orders.iter().sum::<f64>() + f64::from(bracket_hydrogens);

    let eff = if is_constrained(z) {
        effective_atomic_num(z, charge)
    } else {
        z
    };
    let dv = default_valence(eff);

    // An aromatic atom whose 1.5-per-bond sum overshoots is pulled back to the largest valence it
    // does not exceed, provided the overshoot is within one aromatic bond.
    if accum > f64::from(dv) && aromatic {
        let mut pval = f64::from(dv);
        if let Some(valens) = valence::valence_list(eff) {
            for &val in valens {
                if val == -1 || f64::from(val) > accum {
                    break;
                }
                pval = f64::from(val);
            }
        }
        if accum - pval <= AROMATIC_CONTRIB {
            accum = pval;
        }
    }

    accum += 0.1;
    accum.round() as i32
}

/// Implicit hydrogens on an atom, given its explicit valence.
pub fn implicit_valence(
    z: u8,
    charge: i8,
    radicals: u8,
    explicit_valence: i32,
    aromatic: bool,
) -> i32 {
    if z == 0 {
        return 0;
    }
    // a bare hydrogen carries one implicit partner unless charged
    if explicit_valence == 0 && radicals == 0 && z == 1 {
        return match charge {
            1 | -1 => 0,
            0 => 1,
            _ => 0,
        };
    }

    let mut explicit_plus_rad = explicit_valence + i32::from(radicals);

    let mut eff = if is_constrained(z) {
        effective_atomic_num(z, charge)
    } else {
        z
    };
    if eff == 0 {
        return 0;
    }
    let dv = default_valence(eff);
    if dv == -1 {
        return 0;
    }

    if can_be_hypervalent(z, eff) {
        eff = z;
        explicit_plus_rad -= i32::from(charge);
    }
    let Some(valens) = valence::valence_list(eff) else {
        return 0;
    };

    if aromatic {
        // an aromatic atom takes only what its default valence leaves spare
        if explicit_plus_rad <= dv {
            dv - explicit_plus_rad
        } else {
            0
        }
    } else {
        for &v in valens {
            if v < 0 {
                break;
            }
            if explicit_plus_rad <= i32::from(v) {
                return i32::from(v) - explicit_plus_rad;
            }
        }
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neutral_organics() {
        // methane carbon: no bonds, no bracket H -> four implicit
        assert_eq!(implicit_valence(6, 0, 0, 0, false), 4);
        // an amine N with two single bonds -> one implicit
        assert_eq!(implicit_valence(7, 0, 0, 2, false), 1);
        // ether O with two single bonds -> none
        assert_eq!(implicit_valence(8, 0, 0, 2, false), 0);
    }

    #[test]
    fn zero_h_bracket_atoms_are_no_implicit() {
        // `[O]` is atomic oxygen (0 H, a radical), NOT organic `O` (water): brackets mean the H count
        // is stated, so an absent count is zero, not implicit. This exercises the parser fix through
        // counts(): 0 hydrogens and 2 radical electrons, matching RDKit.
        let c = counts(&crate::parse("[O]").unwrap());
        assert_eq!(c.total_num_hs, vec![0]);
        assert_eq!(c.radicals, vec![2]);
        // organic `O` is unchanged — water, 2 implicit H, no radical.
        let w = counts(&crate::parse("O").unwrap());
        assert_eq!(w.total_num_hs, vec![2]);
        assert_eq!(w.radicals, vec![0]);
        // a satisfied charged bracket atom keeps 0 H and 0 radicals (the corpus case, gate-preserving).
        let a = counts(&crate::parse("CC(=O)[O-]").unwrap());
        assert_eq!(a.total_num_hs, vec![3, 0, 0, 0]);
        assert_eq!(a.radicals, vec![0, 0, 0, 0]);
    }

    #[test]
    fn radicals_match_rdkit_assignradicals() {
        // Genuine radicals (bracket atoms whose stated valence leaves unpaired electrons):
        // carbene [CH2] — 0 bonds + 2 H → totalValence 2 → 2 radical electrons (triplet).
        assert_eq!(radical_electrons(6, 0, true, 2, 0), 2);
        // atomic oxygen [O] — 0 bonds, 0 H → 2 radical electrons (triplet O).
        assert_eq!(radical_electrons(8, 0, true, 0, 0), 2);
        // methyl radical [CH3] — 3 H → totalValence 3 → 1 radical electron.
        assert_eq!(radical_electrons(6, 0, true, 3, 0), 1);

        // Closed-shell atoms must stay at zero (this is what keeps the perception gate intact):
        // a non-bracket atom is never assigned radicals (it carries implicit H instead).
        assert_eq!(radical_electrons(6, 0, false, 2, 0), 0);
        // 4-valent bracket stereocentre [C@@H](*)(*)* — 3 bonds + 1 H.
        assert_eq!(radical_electrons(6, 0, true, 4, 3), 0);
        // carboxylate [O-] — 1 bond, charge −1.
        assert_eq!(radical_electrons(8, -1, true, 1, 1), 0);
        // quaternary [N+] — 4 bonds, charge +1.
        assert_eq!(radical_electrons(7, 1, true, 4, 4), 0);
        // pyrrole [nH] — two aromatic bonds (3.0) + 1 H → int(4.1) = 4.
        assert_eq!(radical_electrons(7, 0, true, 4, 2), 0);
    }

    #[test]
    fn charge_shifts_the_element() {
        // quaternary N+ reads carbon's list, so four bonds satisfy it
        assert_eq!(implicit_valence(7, 1, 0, 4, false), 0);
        // carboxylate O- reads nitrogen's list
        assert_eq!(implicit_valence(8, -1, 0, 1, false), 0);
    }

    #[test]
    fn sulfur_reaches_six() {
        // sulfonamide S: two double bonds to O plus two single bonds
        assert_eq!(implicit_valence(16, 0, 0, 6, false), 0);
        // thioether S with two single bonds
        assert_eq!(implicit_valence(16, 0, 0, 2, false), 0);
        // sulfoxide S at 4
        assert_eq!(implicit_valence(16, 0, 0, 4, false), 0);
    }

    #[test]
    fn aromatic_bonds_round() {
        // aromatic CH in benzene: two aromatic bonds -> 3.0, +0.1, rounds to 3
        assert_eq!(explicit_valence(6, 0, 0, true, &[1.5, 1.5]), 3);
        assert_eq!(implicit_valence(6, 0, 0, 3, true), 1);
        // Pyrrole's N: two aromatic bonds (3.0) plus its bracket H is 4.0, which overshoots the
        // default valence of 3 by one aromatic bond, so it is pulled back. RDKit reports explicit
        // valence 3 and implicit 0 here; the bracket H is the one hydrogen.
        assert_eq!(explicit_valence(7, 0, 1, true, &[1.5, 1.5]), 3);
        assert_eq!(implicit_valence(7, 0, 0, 3, true), 0);
        // pyridine's N has no hydrogen
        assert_eq!(explicit_valence(7, 0, 0, true, &[1.5, 1.5]), 3);
        assert_eq!(implicit_valence(7, 0, 0, 3, true), 0);
    }
}

/// Per-atom valence and hydrogen counts for a parsed molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Counts {
    pub explicit_valence: Vec<i32>,
    pub implicit_valence: Vec<i32>,
    pub total_num_hs: Vec<i32>,
    /// Radical electrons per atom (RDKit `assignRadicals`). Zero for atoms that carry implicit
    /// hydrogens (only `noImplicit` / bracket atoms are assigned radicals) and for closed-shell
    /// bracket atoms; nonzero only for genuine radicals (e.g. `[CH2]`, `[O]`).
    pub radicals: Vec<u8>,
}

/// Radical electrons on one atom — a faithful port of RDKit `MolOps::assignRadicals`
/// (`GraphMol/MolOps.cpp`). Only `noImplicit` atoms (those written in brackets, which state their
/// hydrogens) are assigned radicals; every other atom keeps implicit hydrogens and so has none.
/// `total_valence` is the raw bond-order sum plus explicit (bracket) hydrogens, `int(accum + 0.1)`.
fn radical_electrons(z: u8, charge: i8, no_implicit: bool, total_valence: i32, degree: usize) -> u8 {
    if !no_implicit || z == 0 {
        return 0;
    }
    let chg = i32::from(charge);
    let nouter = valence::n_outer_elecs(z).map_or(0, i32::from);
    let valens = valence::valence_list(z);
    // "has a defined valence list" — false only for the `[-1]` sentinel (e.g. transition metals).
    let has_valence_info = valens.is_some_and(|v| !(v.len() == 1 && v[0] == -1));
    if has_valence_info {
        let base_count = if z == 1 || z == 2 { 2 } else { 8 };
        let mut num_radicals = base_count - nouter - total_valence + chg;
        if num_radicals < 0 {
            num_radicals = 0;
            // hypervalent atoms (RDKit github #447): the smallest listed valence that isn't exceeded
            if let Some(v) = valens {
                if v.len() > 1 {
                    for &val in v {
                        if i32::from(val) - total_valence + chg >= 0 {
                            num_radicals = i32::from(val) - total_valence + chg;
                            break;
                        }
                    }
                }
            }
        }
        let num_radicals2 = nouter - total_valence - chg;
        if num_radicals2 >= 0 {
            num_radicals = num_radicals.min(num_radicals2);
        }
        num_radicals.max(0) as u8
    } else {
        // no preferred valence (e.g. transition metals): a lone, bonded metal gets none; an isolated
        // atom carries `(nouter − charge) mod 2` unpaired electrons.
        if degree > 0 {
            0
        } else {
            let n_valence = (nouter - chg).max(0);
            (n_valence % 2) as u8
        }
    }
}

/// Valence contribution of a bond, given whether both ends are aromatic as written.
pub fn contribution(kind: yowl::feature::BondKind, both_aromatic: bool) -> f64 {
    use yowl::feature::BondKind as B;
    match kind {
        B::Aromatic => AROMATIC_CONTRIB,
        B::Double => 2.0,
        B::Triple => 3.0,
        B::Quadruple => 4.0,
        // an unwritten bond between two aromatic atoms is aromatic
        B::Elided if both_aromatic => AROMATIC_CONTRIB,
        _ => 1.0,
    }
}

/// Compute valences and hydrogen counts for every atom of a parsed graph.
/// The per-bond valence contributions and per-atom aromaticity as the SMILES wrote them — RDKit's
/// state at its first `updatePropertyCache` (before kekulization/aromaticity perception).
pub fn as_written_bond_contribs(g: &crate::SmilesGraph) -> Vec<f64> {
    g.bonds
        .iter()
        .enumerate()
        .map(|(bi, &(a, b))| {
            let both = g.atoms[a].aromatic_as_written && g.atoms[b].aromatic_as_written;
            contribution(g.bond_kinds[bi], both)
        })
        .collect()
}

/// Per-atom valence and hydrogen counts, computed from the AS-WRITTEN bond state (RDKit's first
/// `updatePropertyCache`, which feeds kekulization). For the final counts against perceived
/// aromaticity, use [`counts_with`].
pub fn counts(g: &crate::SmilesGraph) -> Counts {
    let contribs = as_written_bond_contribs(g);
    let arom: Vec<bool> = g.atoms.iter().map(|a| a.aromatic_as_written).collect();
    counts_with(g, &contribs, &arom)
}

/// Per-atom valence and hydrogen counts against a given bond state — `bond_contrib[bi]` is bond `bi`'s
/// valence contribution and `atom_aromatic[i]` atom `i`'s aromaticity. RDKit computes this twice: once
/// on the as-written molecule (see [`counts`]) and once after `setAromaticity` on the perceived one;
/// the second is what downstream perception uses. `contribution` per-bond is factored out here so the
/// same routine serves both, mirroring RDKit rewriting bond types before its second pass.
pub fn counts_with(g: &crate::SmilesGraph, bond_contrib: &[f64], atom_aromatic: &[bool]) -> Counts {
    let n = g.atoms.len();
    let mut orders: Vec<Vec<f64>> = vec![Vec::new(); n];
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        // A dative bond contributes 0 to its donor's valence (RDKit's `Bond::DATIVE`); the bond is
        // still present, so it still counts toward degree (`orders[i].len()`).
        orders[a].push(if g.dative_donor[bi] == Some(a) { 0.0 } else { bond_contrib[bi] });
        orders[b].push(if g.dative_donor[bi] == Some(b) { 0.0 } else { bond_contrib[bi] });
    }

    let mut ev_out = Vec::with_capacity(n);
    let mut iv_out = Vec::with_capacity(n);
    let mut hs_out = Vec::with_capacity(n);
    let mut rad_out = Vec::with_capacity(n);
    for (i, atom) in g.atoms.iter().enumerate() {
        let bracket = atom.bracket_hydrogens.unwrap_or(0);
        let no_implicit = atom.bracket_hydrogens.is_some();
        let ev = explicit_valence(
            atom.atomic_number,
            atom.charge,
            bracket,
            atom_aromatic[i],
            &orders[i],
        );
        // a bracket atom states its hydrogens, so none are implied
        let iv = if no_implicit {
            0
        } else {
            implicit_valence(
                atom.atomic_number,
                atom.charge,
                0,
                ev,
                atom_aromatic[i],
            )
        };
        // RDKit's assignRadicals totalValence: the raw bond-order sum plus explicit (bracket) Hs,
        // `int(accum + 0.1)` — NOT the aromatic-adjusted `explicit_valence`.
        let raw_total_valence = (orders[i].iter().sum::<f64>() + f64::from(bracket) + 0.1) as i32;
        rad_out.push(radical_electrons(
            atom.atomic_number,
            atom.charge,
            no_implicit,
            raw_total_valence,
            orders[i].len(),
        ));
        ev_out.push(ev);
        iv_out.push(iv);
        hs_out.push(iv + i32::from(bracket));
    }
    Counts {
        explicit_valence: ev_out,
        implicit_valence: iv_out,
        total_num_hs: hs_out,
        radicals: rad_out,
    }
}
