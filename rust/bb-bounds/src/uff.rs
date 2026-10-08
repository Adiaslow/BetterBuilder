//! UFF atom typing: the per-atom type label and its force-field parameters.
//!
//! A faithful port of RDKit's `RDKit::UFF::Tools::getAtomLabel` + `addAtomChargeFlags`
//! (`Code/GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp`) and the parameter table
//! (`ForceField/UFF/Params.cpp`). `setTopolBounds` reads these params — `r1` gives 1-2 bond rest
//! lengths, `GMP_Xi`/`Z1` the electronegativity correction — so the typing is the correctness
//! foundation the whole bounds matrix sits on.
//!
//! The label is built from element symbol + a hybridization suffix (`R` for conjugated C/N/O/S)
//! and an oxidation-state suffix, then looked up in the parameter table. Every element branch of
//! the source is reproduced, not only those the current corpus exercises. The one RDKit case with
//! no representation here is `SP2D` (square-planar, label suffix `4`): `bb_perceive`'s
//! [`Hybridization`] does not model it, so it cannot arise; it is called out below rather than
//! silently mapped to something else.

use bb_perceive::hybrid::Hybridization;
use bb_perceive::valence;
use std::collections::HashMap;
use std::sync::OnceLock;

/// One row of the UFF parameter table. Field names follow `ForceFields::UFF::AtomicParams`.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomicParams {
    pub r1: f64,
    pub theta0: f64,
    pub x1: f64,
    pub d1: f64,
    pub zeta: f64,
    pub z1: f64,
    pub v1: f64,
    pub u1: f64,
    pub gmp_xi: f64,
    pub gmp_hardness: f64,
    pub gmp_radius: f64,
}

/// RDKit's UFF parameter table, vendored verbatim from `ForceFields::UFF::defaultParamData`.
const PARAM_TABLE: &str = include_str!("../../../rdkit-patch/uff_params.txt");

fn params() -> &'static HashMap<String, AtomicParams> {
    static TABLE: OnceLock<HashMap<String, AtomicParams>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut m = HashMap::new();
        for line in PARAM_TABLE.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            // whitespace-separated, matching the boost tokenizer RDKit parses the table with
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 12 {
                continue;
            }
            let g = |i: usize| f[i].parse::<f64>().unwrap_or(0.0);
            m.insert(
                f[0].to_string(),
                AtomicParams {
                    r1: g(1),
                    theta0: g(2),
                    x1: g(3),
                    d1: g(4),
                    zeta: g(5),
                    z1: g(6),
                    v1: g(7),
                    u1: g(8),
                    gmp_xi: g(9),
                    gmp_hardness: g(10),
                    gmp_radius: g(11),
                },
            );
        }
        m
    })
}

/// The parameters for a UFF type label, or `None` if the table has no such type (RDKit's
/// `foundAll = false`).
pub fn params_for_label(label: &str) -> Option<&'static AtomicParams> {
    params().get(label)
}

/// `Params::lambda`, the scaling factor for the Pauling bond-order correction.
const LAMBDA: f64 = 0.1332;

/// The UFF bond rest length between two atoms, reproducing
/// `ForceFields::UFF::Utils::calcBondRestLength`.
///
/// `bond_order` is `getBondTypeAsDouble` (1.5 for aromatic). The result is the sum of the two bond
/// radii, a Pauling bond-order shortening (`rBO`), and an O'Keeffe–Brese electronegativity
/// correction (`rEN`): `res = ri + rj + rBO - rEN`.
pub fn calc_bond_rest_length(bond_order: f64, end1: &AtomicParams, end2: &AtomicParams) -> f64 {
    let (ri, rj) = (end1.r1, end2.r1);
    // Pauling bond-order correction
    let r_bo = -LAMBDA * (ri + rj) * bond_order.ln();
    // O'Keeffe–Brese electronegativity correction
    let (xi, xj) = (end1.gmp_xi, end2.gmp_xi);
    let d = xi.sqrt() - xj.sqrt();
    let r_en = ri * rj * d * d / (xi * ri + xj * rj);
    ri + rj + r_bo - r_en
}

/// RDKit's `PeriodicTable::getDefaultValence` = the first entry of the valence list (`-1` if none).
fn default_valence(z: u8) -> i32 {
    valence::element(z)
        .and_then(|e| e.valences.first())
        .map_or(-1, |&v| i32::from(v))
}

/// The hybridization suffix for the general (non-special-element) branch, or `None` for the
/// hybridizations RDKit's switch leaves without a suffix (`S`, and the unspecified default).
fn hybridization_suffix(
    hyb: Hybridization,
    z: u8,
    aromatic: bool,
    has_conjugated_bond: bool,
) -> Option<char> {
    match hyb {
        // RDKit `case S:` does nothing; UNSPECIFIED falls to the default (error, no suffix).
        Hybridization::Unspecified | Hybridization::S => None,
        Hybridization::Sp => Some('1'),
        Hybridization::Sp2 => {
            // aromatic or conjugated C/N/O/S becomes the resonant type `R`
            if (aromatic || has_conjugated_bond) && matches!(z, 6 | 7 | 8 | 16) {
                Some('R')
            } else {
                Some('2')
            }
        }
        Hybridization::Sp3 => Some('3'),
        // RDKit `SP2D` -> '4' has no bb_perceive Hybridization and cannot arise here.
        Hybridization::Sp3d => Some('5'),  // RDKit SP3D
        Hybridization::Sp3d2 => Some('6'), // RDKit SP3D2
    }
}

/// Append RDKit's oxidation-state charge flag to `key`, reproducing `addAtomChargeFlags`.
///
/// `tolerate_charge_mismatch` is `true` at the `getAtomLabel` call site (the header default), so a
/// valence that does not match the expected state still gets the fallback suffix. Elements not
/// listed here (H, C, N, O, halogens, …) get no suffix. Al and Si get none by design.
fn add_charge_flags(
    z: u8,
    total_valence: i32,
    formal_charge: i32,
    hyb: Hybridization,
    key: &mut String,
) {
    let tolerate = true;
    let tv = total_valence;
    let fc = formal_charge;
    // append `suf` when the valence matches `want` or `fc` matches, else when tolerating
    let single = |want: i32, suf: &str, key: &mut String| {
        if tv == want || fc == want || tolerate {
            key.push_str(suf);
        }
    };
    match z {
        29 | 47 => single(1, "+1", key),                         // Cu, Ag
        4 | 20 | 25 | 26 | 28 | 46 | 78 => single(2, "+2", key), // Be, Ca, Mn, Fe, Ni, Pd, Pt
        21 | 24 | 27 | 79 | 89 | 96 | 97 | 98 | 99 | 100 | 101 | 102 | 103 => single(3, "+3", key),
        2 | 18 | 22 | 36 | 54 | 90 | 91 | 92 | 93 | 94 | 95 => single(4, "+4", key),
        23 | 41 | 43 | 73 => single(5, "+5", key), // V, Nb, Tc, Ta
        42 => single(6, "+6", key),                // Mo
        12 => key.push_str("+2"),                  // Mg: 2 -> +2, else (tolerate) +2
        13 => {}                                   // Al: warn only, no suffix
        14 => {}                                   // Si: warn only, no suffix
        15 => match tv {
            // P
            3 => key.push_str("+3"),
            5 => key.push_str("+5"),
            _ => key.push_str("+5"), // tolerate fallback
        },
        // S: only when not SP2
        16 if hyb != Hybridization::Sp2 => match tv {
            2 => key.push_str("+2"),
            4 => key.push_str("+4"),
            6 => key.push_str("+6"),
            _ => key.push_str("+6"), // tolerate fallback
        },
        30 | 48 | 34 | 52 | 80 | 84 => key.push_str("+2"), // Zn, Cd, Se, Te, Hg, Po
        31 | 33 | 49 | 51 | 81 | 82 | 83 => key.push_str("+3"), // Ga, As, In, Sb, Tl, Pb, Bi
        75 => {
            // Re: only via the tolerate path, and it rewrites the key
            if key == "Re6" {
                *key = "Re6+5".to_string();
            } else if key == "Re3" {
                *key = "Re3+7".to_string();
            }
        }
        _ => {}
    }
    // lanthanides: 6 -> +3, else (tolerate) +3
    if (57..=71).contains(&z) {
        key.push_str("+3");
    }
}

/// The UFF atom-type label for one atom, reproducing `getAtomLabel`.
///
/// `total_valence` is `getTotalValence`, `has_conjugated_bond` is `atomHasConjugatedBond`; both are
/// carried on [`bb_perceive::smarts_match::Perceived`].
pub fn atom_label(
    z: u8,
    hyb: Hybridization,
    aromatic: bool,
    has_conjugated_bond: bool,
    total_valence: i32,
    formal_charge: i32,
) -> String {
    let symbol = valence::element(z).map_or("*", |e| e.symbol.as_str());
    let mut key = symbol.to_string();
    if key.len() == 1 {
        key.push('_');
    }

    if z != 0 {
        let nouter = valence::element(z).map_or(0, |e| i32::from(e.n_outer_elecs));
        // hybridization is skipped for alkali metals (1 outer e-) and halogens (7 outer e-), unless
        // the element has no default valence at all
        if default_valence(z) == -1 || (nouter != 1 && nouter != 7) {
            match z {
                // main-group elements UFF always types as their '3' / '1' state
                12 | 13 | 14 | 15 | 50 | 51 | 52 | 81 | 82 | 83 | 84 => key.push('3'),
                80 => key.push('1'), // Hg
                _ => {
                    if let Some(s) = hybridization_suffix(hyb, z, aromatic, has_conjugated_bond) {
                        key.push(s);
                    }
                }
            }
        }
    }

    add_charge_flags(z, total_valence, formal_charge, hyb, &mut key);
    key
}

/// One bond's rest length together with whether both its atoms had UFF params.
///
/// `set12Bounds` sets *different* distance bounds in the two cases (`bl ± tol` when params exist,
/// `1.5·bl`/`0.5·bl` in the fallback), so the branch has to be visible, not just the length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondLen {
    pub length: f64,
    pub has_params: bool,
}

/// Per-bond UFF rest length and params-present flag, in bond order — `set12Bounds`'
/// `accumData.bondLengths` plus the branch it took.
///
/// The length is `calc_bond_rest_length` from the two atoms' UFF params when both have params and
/// the bond order is positive, and RDKit's crude `(rvdw1 + rvdw2) / 2` fallback otherwise.
pub fn bond_rest_lengths_detailed(mol: &bb_perceive::smarts_match::Perceived) -> Vec<BondLen> {
    let labels = atom_labels(mol);
    mol.bonds
        .iter()
        .enumerate()
        .map(|(bi, &(a, b))| {
            // getBondTypeAsDouble: 1.5 for aromatic, otherwise the integer order
            let bo = if mol.bond_aromatic[bi] {
                1.5
            } else {
                f64::from(mol.bond_order[bi])
            };
            match (params_for_label(&labels[a]), params_for_label(&labels[b])) {
                (Some(pa), Some(pb)) if bo > 0.0 => BondLen {
                    length: calc_bond_rest_length(bo, pa, pb),
                    has_params: true,
                },
                _ => {
                    // no UFF params for an atom: RDKit falls back to the mean vdW radius
                    let vw1 = valence::rvdw(mol.atomic_numbers[a]).unwrap_or(0.0);
                    let vw2 = valence::rvdw(mol.atomic_numbers[b]).unwrap_or(0.0);
                    BondLen {
                        length: (vw1 + vw2) / 2.0,
                        has_params: false,
                    }
                }
            }
        })
        .collect()
}

/// Per-bond UFF rest length only — `set12Bounds`' `accumData.bondLengths`.
pub fn bond_rest_lengths(mol: &bb_perceive::smarts_match::Perceived) -> Vec<f64> {
    bond_rest_lengths_detailed(mol)
        .iter()
        .map(|b| b.length)
        .collect()
}

/// UFF out-of-plane (inversion) coefficients and force constant for an sp2 center, reproducing
/// `ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant`. Returns `(res, c0, c1, c2)`
/// where `res` is the pre-scaling force constant.
pub fn calc_inversion_coefficients(center_z: u8, is_c_bound_to_o: bool) -> (f64, f64, f64, f64) {
    let (mut res, c0, c1, c2);
    if center_z == 6 || center_z == 7 || center_z == 8 {
        c0 = 1.0;
        c1 = -1.0;
        c2 = 0.0;
        res = if is_c_bound_to_o { 50.0 } else { 6.0 };
    } else {
        // group-5 centers (P/As/Sb/Bi): the UFF paper is unclear; these constants follow RDKit
        let mut w0 = std::f64::consts::PI / 180.0;
        w0 *= match center_z {
            15 => 84.4339,
            33 => 86.9735,
            51 => 87.7047,
            83 => 90.0,
            _ => 0.0,
        };
        c2 = 1.0;
        c1 = -4.0 * w0.cos();
        c0 = -(c1 * w0.cos() + c2 * (2.0 * w0).cos());
        res = 22.0 / (c0 + c1 + c2);
    }
    res /= 3.0;
    (res, c0, c1, c2)
}

/// Per-atom UFF type labels for a perceived molecule, in its (post-`AddHs`) atom order.
pub fn atom_labels(mol: &bb_perceive::smarts_match::Perceived) -> Vec<String> {
    (0..mol.atomic_numbers.len())
        .map(|i| {
            atom_label(
                mol.atomic_numbers[i],
                mol.hybridization[i],
                mol.aromatic[i],
                mol.atom_conjugated[i],
                mol.total_valence[i],
                i32::from(mol.charges[i]),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn organic_labels() {
        // sp3 carbon -> C_3; aromatic carbon -> C_R; carbonyl carbon (sp2, conjugated) -> C_R
        assert_eq!(atom_label(6, Hybridization::Sp3, false, false, 4, 0), "C_3");
        assert_eq!(atom_label(6, Hybridization::Sp2, true, true, 4, 0), "C_R");
        // non-conjugated sp2 carbon (e.g. an isolated alkene) -> C_2
        assert_eq!(atom_label(6, Hybridization::Sp2, false, false, 4, 0), "C_2");
        // hydrogen: no hybridization suffix, no charge flag
        assert_eq!(
            atom_label(1, Hybridization::Unspecified, false, false, 1, 0),
            "H_"
        );
        // amide nitrogen (sp2, conjugated) -> N_R; amine nitrogen (sp3) -> N_3
        assert_eq!(atom_label(7, Hybridization::Sp2, false, true, 3, 0), "N_R");
        assert_eq!(atom_label(7, Hybridization::Sp3, false, false, 3, 0), "N_3");
        // carbonyl / hydroxyl oxygen
        assert_eq!(atom_label(8, Hybridization::Sp2, false, true, 2, 0), "O_R");
        assert_eq!(atom_label(8, Hybridization::Sp3, false, false, 2, 0), "O_3");
    }

    #[test]
    fn halogens_skip_hybridization() {
        assert_eq!(atom_label(9, Hybridization::Sp3, false, false, 1, 0), "F_");
        assert_eq!(atom_label(17, Hybridization::Sp3, false, false, 1, 0), "Cl");
        assert_eq!(atom_label(35, Hybridization::Sp3, false, false, 1, 0), "Br");
    }

    #[test]
    fn sulfur_and_phosphorus_oxidation_states() {
        // sulfone / sulfonamide S (valence 6, not sp2) -> S_3+6
        assert_eq!(
            atom_label(16, Hybridization::Sp3, false, false, 6, 0),
            "S_3+6"
        );
        // thioether S (valence 2) -> S_3+2
        assert_eq!(
            atom_label(16, Hybridization::Sp3, false, false, 2, 0),
            "S_3+2"
        );
        // aromatic/sp2 S takes no charge flag -> S_R
        assert_eq!(atom_label(16, Hybridization::Sp2, true, true, 2, 0), "S_R");
        // phosphate P (valence 5) -> P_3+5; phosphine P (valence 3) -> P_3+3
        assert_eq!(
            atom_label(15, Hybridization::Sp3, false, false, 5, 0),
            "P_3+5"
        );
        assert_eq!(
            atom_label(15, Hybridization::Sp3, false, false, 3, 0),
            "P_3+3"
        );
    }

    #[test]
    fn bond_rest_length_matches_hand_computation() {
        // A C_3-C_3 single bond: ln(1)=0 kills rBO, and equal electronegativities kill rEN, so the
        // rest length is just 2*r1(C_3) = 2*0.757 = 1.514.
        let c = params_for_label("C_3").unwrap();
        assert!((calc_bond_rest_length(1.0, c, c) - 1.514).abs() < 1e-12);
        // A double bond between the same atoms shortens via rBO = -lambda*(ri+rj)*ln(2).
        let expected = 1.514 - 0.1332 * (0.757 + 0.757) * 2.0f64.ln();
        assert!((calc_bond_rest_length(2.0, c, c) - expected).abs() < 1e-12);
    }

    #[test]
    fn missing_params_use_vdw_fallback() {
        use bb_perceive::smarts_match::Perceived;
        // a dummy atom (atomic number 0, UFF label "*", no params) single-bonded to an sp3 carbon:
        // set12's else-branch applies, (rvdw(0) + rvdw(6)) / 2
        let mol = Perceived {
            atomic_numbers: vec![0, 6],
            charges: vec![0, 0],
            aromatic: vec![false, false],
            hybridization: vec![Hybridization::Unspecified, Hybridization::Sp3],
            bonds: vec![(0, 1)],
            bond_order: vec![1],
            bond_aromatic: vec![false],
            bond_in_ring: vec![false],
            rings: vec![],
            total_valence: vec![0, 4],
            atom_conjugated: vec![false, false],
            bond_conjugated: vec![false],
            bond_stereo: vec![bb_perceive::smarts_match::Stereo::None],
            stereo_atoms: vec![[-1, -1]],
            chiral_tags: vec![bb_perceive::chirality::ChiralTag::None; 2],
        };
        assert!(
            params_for_label("*").is_none(),
            "dummy atom must have no UFF params"
        );
        let bl = bond_rest_lengths(&mol);
        let expected = (bb_perceive::valence::rvdw(0).unwrap_or(0.0)
            + bb_perceive::valence::rvdw(6).unwrap())
            / 2.0;
        assert!(
            (bl[0] - expected).abs() < 1e-12,
            "fallback {} != {expected}",
            bl[0]
        );
    }

    #[test]
    fn every_produced_label_resolves_in_the_table() {
        for lbl in [
            "C_3", "C_R", "C_2", "H_", "N_R", "N_3", "O_R", "O_3", "F_", "Cl", "Br", "S_3+6",
            "S_3+2", "S_R", "P_3+5", "P_3+3",
        ] {
            assert!(params_for_label(lbl).is_some(), "table missing {lbl}");
        }
    }
}
