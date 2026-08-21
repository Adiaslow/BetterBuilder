//! bb-rdkit — setup: SMILES → [`bb_core::MoleculeSpec`], via `cxx` FFI to the patched RDKit C++
//! (2026.09.1pre + the amide-torsion patch).
//!
//! [`build_spec`] extracts the bounds matrix, experimental torsions, chiral and tetrahedral sets,
//! UFF impropers, bonds and angles, double-bond stereo, and the core-pin recipe counts.

use bb_core::{Angle, ChiralSet, ExpTorsion, Improper, MoleculeSpec, StereoDoubleBond};
use cxx::vector::VectorElement;
use cxx::CxxVector;
use thiserror::Error;

/// Errors from setup extraction. [`SpecError::MalformedFfi`] is returned when an array from the C++
/// bridge is shorter than its declared per-item stride implies.
#[derive(Debug, Error)]
pub enum SpecError {
    #[error("RDKit could not parse SMILES: {0:?}")]
    InvalidSmiles(String),
    #[error("malformed data from the RDKit FFI bridge: index {0} out of range")]
    MalformedFfi(usize),
    #[error("RDKit threw while perceiving {smiles:?}: {message}")]
    RdkitThrew { smiles: String, message: String },
}

/// Copy element `i` out of an FFI vector, or [`SpecError::MalformedFfi`] if out of range.
#[inline]
fn at<T: VectorElement + Copy>(v: &CxxVector<T>, i: usize) -> Result<T, SpecError> {
    v.get(i).copied().ok_or(SpecError::MalformedFfi(i))
}

#[cxx::bridge(namespace = "bb")]
mod ffi {
    unsafe extern "C++" {
        include!("bb-rdkit/src/bridge.h");

        fn num_heavy_atoms(smiles: &str) -> Result<usize>;

        /// RDKit's `triangleSmoothBounds` on a flat n×n bounds matrix.
        fn smooth_bounds(bounds: &[f64], n: usize, tol: f64) -> Result<Vec<f64>>;

        /// RDKit's raw `setTopolBounds` matrix (pre-smoothing), flat n×n. Empty on parse failure.
        fn raw_bounds(smiles: &str) -> Result<Vec<f64>>;

        /// RDKit's UFF parameter table (`defaultParamData`) verbatim, for vendoring.
        fn uff_param_table() -> Result<String>;

        /// Per-atom UFF atom-type labels in post-AddHs order. Empty on parse failure.
        fn uff_atom_labels(smiles: &str) -> Result<Vec<String>>;

        /// Per-bond UFF rest lengths (set12's accumData.bondLengths), post-AddHs bond order.
        fn uff_bond_rest_lengths(smiles: &str) -> Result<Vec<f64>>;

        /// RDKit's topological distance matrix (getDistanceMat), flat n×n. Empty on parse failure.
        fn topo_distance_matrix(smiles: &str) -> Result<Vec<f64>>;

        /// RDKit's PERCEIVED per-atom + per-bond aromaticity (post-sanitize, heavy-atom/SMILES order):
        /// [n_atoms, atomArom×n_atoms, (begin,end,bondArom)×n_bonds]; leading -1 = parse failure.
        fn set_arom_only(smiles: &str) -> Result<Vec<i32>>;
        fn atom_elec(smiles: &str) -> Result<Vec<i32>>;
        fn aromatic_perception(smiles: &str) -> Result<Vec<i32>>;

        /// RDKit's per-bond double-bond stereo: 3 ints per bond (stereo enum, 2 stereo atoms).
        fn bond_dirs(smiles: &str) -> Result<Vec<i32>>;
        /// RDKit Stage-A DistGeom force-field energy at 4D coords (firstMinimization FF). NaN on failure.
        fn stage_a_ff_grad(smiles: &str, coords: &[f64], weight_chiral: f64, weight_fourth: f64, basin: f64) -> Vec<f64>;
        fn stage_c_ff_grad(smiles: &str, coords: &[f64]) -> Vec<f64>;
        fn stage_c_improper_grad(smiles: &str, coords: &[f64]) -> Vec<f64>;
        fn stage_c_torsion_grad(smiles: &str, coords: &[f64]) -> Vec<f64>;
        /// RDKit's post-embed accept/reject checks on a 3D conformer, in native passes_checks order:
        /// 6 bits {tetrahedral, chiral, planarity, double-bond-geometry, final-chiral,
        /// double-bond-stereo}, 1=pass. Empty on failure.
        fn embed_checks(smiles: &str, coords: &[f64]) -> Vec<i32>;
        /// RDKit's coordMap-tightened bounds (setupInitialBoundsMatrix success path), n*n row-major.
        /// `pin_idx[k]` = atom index of pin k; `pin_xyz` = 3 f64 per pin. Empty on failure.
        fn coord_map_bounds(smiles: &str, pin_idx: &[i32], pin_xyz: &[f64]) -> Vec<f64>;
        /// RDKit's powerEigenSolver on a metric matrix T (lower triangle incl. diagonal, size
        /// n(n+1)/2). Returns [n_eig eigenvalues][n_eig*n eigvec components]. Empty on size mismatch.
        /// RDKit's Stage-A calcEnergy at a 4D geometry (the value the per-atom reject thresholds).
        fn stage_a_energy(smiles: &str, coords: &[f64], weight_chiral: f64, weight_fourth: f64, basin: f64) -> f64;
        fn bond_types(smiles: &str) -> Result<Vec<i32>>;
        fn sssr_rings(smiles: &str) -> Result<Vec<i32>>;
        fn bond_stereo(smiles: &str) -> Result<Vec<i32>>;

        /// RDKit's per-atom getChiralTag (0=unspecified, 1=CW, 2=CCW), post-AddHs order.
        fn chiral_tags(smiles: &str) -> Result<Vec<i32>>;

        /// RDKit's raw bounds with set13/set14/set15 selectable (false,false,false = set12+VDW).
        fn raw_bounds_stage(smiles: &str, set13: bool, set14: bool, set15: bool) -> Result<Vec<f64>>;

        /// Opaque holder for one molecule's extracted distance-geometry problem.
        type MolSpecC;
        /// `Result` so that cxx wraps the call in a try/catch: RDKit throws on some inputs, and an
        /// exception crossing an unwrapped boundary calls `std::terminate`.
        fn build_mol_spec(smiles: &str) -> Result<UniquePtr<MolSpecC>>;
        fn ok(self: &MolSpecC) -> bool;
        fn n_atoms(self: &MolSpecC) -> usize;
        fn atomic_nums(self: &MolSpecC) -> &CxxVector<u8>; // 1 per atom, addHs order
        fn formal_charge(self: &MolSpecC) -> i32; // net formal charge
        fn bounds(self: &MolSpecC) -> &CxxVector<f64>; // n*n row-major, upper=UB lower=LB
        fn tors_atoms(self: &MolSpecC) -> &CxxVector<u32>; // 4 per torsion
        fn tors_v(self: &MolSpecC) -> &CxxVector<f64>; // 6 per torsion
        fn tors_signs(self: &MolSpecC) -> &CxxVector<i32>; // 6 per torsion
        fn chi_center(self: &MolSpecC) -> &CxxVector<u32>; // 1 per chiral set
        fn chi_atoms(self: &MolSpecC) -> &CxxVector<u32>; // 4 per set
        fn chi_vol(self: &MolSpecC) -> &CxxVector<f64>; // 2 per set (lo, hi)
        fn chi_fused(self: &MolSpecC) -> &CxxVector<u8>; // 1 per set
        fn imp_atoms(self: &MolSpecC) -> &CxxVector<u32>; // 4 per improper contrib
        fn imp_coef(self: &MolSpecC) -> &CxxVector<f64>; // 4 per contrib (c0, c1, c2, fc)
        fn bonds(self: &MolSpecC) -> &CxxVector<u32>; // 2 per bond (begin, end)
        fn angles(self: &MolSpecC) -> &CxxVector<u32>; // 4 per angle (i, j, k, tripleFlag)
        fn pin_atoms(self: &MolSpecC) -> &CxxVector<u32>; // core-pin recipe frozen atoms
        fn core_seeds(self: &MolSpecC) -> u32; // rdkit_confs_1
        fn sidechain_confs(self: &MolSpecC) -> u32; // per-seed confs
        fn dbe(self: &MolSpecC) -> &CxxVector<u32>; // 3 per double-bond end
        fn sdb_atoms(self: &MolSpecC) -> &CxxVector<u32>; // 4 per stereo double bond
        fn sdb_sign(self: &MolSpecC) -> &CxxVector<i32>; // 1 per stereo double bond
    }
}

/// Heavy-atom count via RDKit's SMILES parser. 0 on parse failure or if RDKit throws.
pub fn num_heavy_atoms(smiles: &str) -> usize {
    ffi::num_heavy_atoms(smiles).unwrap_or(0)
}

/// RDKit's `triangleSmoothBounds` on a flat `n×n` bounds matrix (row<col = upper, row>col = lower).
/// Not used by [`build_spec`].
pub fn smooth_bounds(bounds: &[f64], n: usize, tol: f64) -> Result<Vec<f64>, SpecError> {
    ffi::smooth_bounds(bounds, n, tol).map_err(|e| SpecError::RdkitThrew {
        smiles: String::new(),
        message: e.what().to_string(),
    })
}

/// RDKit's raw `setTopolBounds` matrix (before `triangleSmoothBounds`), flat `n×n` row-major
/// (row<col = upper bound, row>col = lower bound). The oracle for the pure-Rust `setTopolBounds`
/// port: it isolates the constructor, since [`build_spec`] only exposes the smoothed matrix.
/// Returns an empty vector for a SMILES RDKit cannot parse.
pub fn raw_bounds(smiles: &str) -> Result<Vec<f64>, SpecError> {
    ffi::raw_bounds(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// RDKit's UFF parameter table (`ForceFields::UFF::defaultParamData`) verbatim — the exact string
/// the UFF atom typer parses. Vendored so the Rust port reads a byte-identical copy.
pub fn uff_param_table() -> Result<String, SpecError> {
    ffi::uff_param_table().map_err(|e| SpecError::RdkitThrew {
        smiles: String::new(),
        message: e.what().to_string(),
    })
}

/// Per-atom UFF atom-type labels (`getAtomLabel`) in post-AddHs order — the oracle for the Rust
/// atom-typing port. Empty for a SMILES RDKit cannot parse.
pub fn uff_atom_labels(smiles: &str) -> Result<Vec<String>, SpecError> {
    ffi::uff_atom_labels(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// Per-bond UFF rest lengths as `set12Bounds` computes them (`accumData.bondLengths`), in post-AddHs
/// bond order — the oracle for the Rust `calc_bond_rest_length` port. Empty on a parse failure.
pub fn uff_bond_rest_lengths(smiles: &str) -> Result<Vec<f64>, SpecError> {
    ffi::uff_bond_rest_lengths(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// RDKit's topological distance matrix (`getDistanceMat`, flat `n×n` row-major, post-AddHs atoms) —
/// the oracle for the Rust `distance_matrix` port. Empty for a SMILES RDKit cannot parse.
pub fn topo_distance_matrix(smiles: &str) -> Result<Vec<f64>, SpecError> {
    ffi::topo_distance_matrix(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// RDKit's per-bond double-bond stereo (`getStereo` enum + two `getStereoAtoms`), 3 ints per bond in
/// post-AddHs bond order — the oracle for the Rust double-bond-stereo port. Empty on a parse failure.
pub fn stage_a_ff_grad(smiles: &str, coords: &[f64], weight_chiral: f64, weight_fourth: f64, basin: f64) -> Vec<f64> {
    ffi::stage_a_ff_grad(smiles, coords, weight_chiral, weight_fourth, basin)
}
pub fn stage_c_ff_grad(smiles: &str, coords: &[f64]) -> Vec<f64> {
    ffi::stage_c_ff_grad(smiles, coords)
}
pub fn stage_c_improper_grad(smiles: &str, coords: &[f64]) -> Vec<f64> {
    ffi::stage_c_improper_grad(smiles, coords)
}
pub fn stage_c_torsion_grad(smiles: &str, coords: &[f64]) -> Vec<f64> {
    ffi::stage_c_torsion_grad(smiles, coords)
}
/// RDKit's five post-embed accept/reject checks on a 3D conformer (see the FFI decl). Empty on failure.
pub fn embed_checks(smiles: &str, coords: &[f64]) -> Vec<i32> {
    ffi::embed_checks(smiles, coords)
}
/// RDKit's coordMap-tightened bounds matrix (see the FFI decl). Empty on failure.
pub fn coord_map_bounds(smiles: &str, pin_idx: &[i32], pin_xyz: &[f64]) -> Vec<f64> {
    ffi::coord_map_bounds(smiles, pin_idx, pin_xyz)
}
/// RDKit's Stage-A calcEnergy at a 4D geometry (firstMinimization FF). NaN on failure.
pub fn stage_a_energy(smiles: &str, coords: &[f64], weight_chiral: f64, weight_fourth: f64, basin: f64) -> f64 {
    ffi::stage_a_energy(smiles, coords, weight_chiral, weight_fourth, basin)
}
pub fn bond_dirs(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::bond_dirs(smiles).map_err(|e| SpecError::RdkitThrew { smiles: smiles.to_string(), message: e.what().to_string() })
}
pub fn bond_types(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::bond_types(smiles).map_err(|e| SpecError::RdkitThrew { smiles: smiles.to_string(), message: e.what().to_string() })
}
pub fn sssr_rings(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::sssr_rings(smiles).map_err(|e| SpecError::RdkitThrew { smiles: smiles.to_string(), message: e.what().to_string() })
}
pub fn bond_stereo(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::bond_stereo(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// RDKit's PERCEIVED aromaticity (post-sanitize) — the oracle for bb-perceive's aromaticity port.
/// `[n_atoms, atomArom×n_atoms, (begin,end,bondArom)×n_bonds]`; a leading `-1` means RDKit rejected
/// the SMILES.
pub fn set_arom_only(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::set_arom_only(smiles).map_err(|e| SpecError::RdkitThrew { smiles: smiles.to_string(), message: e.what().to_string() })
}
pub fn atom_elec(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::atom_elec(smiles).map_err(|e| SpecError::RdkitThrew { smiles: smiles.to_string(), message: e.what().to_string() })
}

pub fn aromatic_perception(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::aromatic_perception(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// RDKit's per-atom `getChiralTag` (0=unspecified, 1=CW, 2=CCW) in post-AddHs order — the oracle for
/// the Rust tetrahedral-chirality port. Empty for a SMILES RDKit cannot parse.
pub fn chiral_tags(smiles: &str) -> Result<Vec<i32>, SpecError> {
    ffi::chiral_tags(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// RDKit's raw bounds matrix with `set13`/`set14`/`set15` selectable — the staged oracle for the
/// bounds-matrix port. `(false, false, false)` is set12 + VDW only. Empty on a parse failure.
pub fn raw_bounds_stage(
    smiles: &str,
    set13: bool,
    set14: bool,
    set15: bool,
) -> Result<Vec<f64>, SpecError> {
    ffi::raw_bounds_stage(smiles, set13, set14, set15).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })
}

/// Build a [`MoleculeSpec`] from a SMILES. Atom indices are in RDKit `SmilesToMol` + `AddHs` order.
///
/// Returns [`SpecError::InvalidSmiles`] if RDKit cannot parse the SMILES, and
/// [`SpecError::RdkitThrew`] if it throws during perception.
pub fn build_spec(smiles: &str) -> Result<MoleculeSpec, SpecError> {
    let c = ffi::build_mol_spec(smiles).map_err(|e| SpecError::RdkitThrew {
        smiles: smiles.to_string(),
        message: e.what().to_string(),
    })?;
    if !c.ok() {
        return Err(SpecError::InvalidSmiles(smiles.to_string()));
    }
    let n = c.n_atoms();
    let bounds: Vec<f32> = c.bounds().iter().map(|&x| x as f32).collect();

    // per-atom element; one entry per atom, so a short vector is a bridge bug
    let atomic_numbers: Vec<u8> = c.atomic_nums().iter().copied().collect();
    if atomic_numbers.len() != n {
        return Err(SpecError::MalformedFfi(atomic_numbers.len()));
    }

    // experimental torsions
    let (ta, tv, ts) = (c.tors_atoms(), c.tors_v(), c.tors_signs());
    let mut exp_torsions = Vec::with_capacity(ta.len() / 4);
    for i in 0..ta.len() / 4 {
        let atoms = [
            at(ta, 4 * i)?,
            at(ta, 4 * i + 1)?,
            at(ta, 4 * i + 2)?,
            at(ta, 4 * i + 3)?,
        ];
        let mut v = [0f64; 6];
        let mut signs = [0i8; 6];
        for k in 0..6 {
            v[k] = at(tv, 6 * i + k)?;
            signs[k] = at(ts, 6 * i + k)? as i8;
        }
        exp_torsions.push(ExpTorsion { atoms, v, signs });
    }

    // chiral sets — vol (0,0) => tetrahedral (checks only), else energy-bearing chiral center
    let (cc, ca, cv, cf) = (c.chi_center(), c.chi_atoms(), c.chi_vol(), c.chi_fused());
    let mut chiral_sets = Vec::new();
    let mut tetrahedral_centers = Vec::new();
    for i in 0..cc.len() {
        let cs = ChiralSet {
            center: at(cc, i)?,
            atoms: [
                at(ca, 4 * i)?,
                at(ca, 4 * i + 1)?,
                at(ca, 4 * i + 2)?,
                at(ca, 4 * i + 3)?,
            ],
            vol_lo: at(cv, 2 * i)?,
            vol_hi: at(cv, 2 * i + 1)?,
            fused_small_rings: at(cf, i)? != 0,
        };
        if cs.vol_lo == 0.0 && cs.vol_hi == 0.0 {
            tetrahedral_centers.push(cs);
        } else {
            chiral_sets.push(cs);
        }
    }

    // UFF impropers — 4 atoms + (c0, c1, c2, fc) per contrib
    let (ia, ic) = (c.imp_atoms(), c.imp_coef());
    let mut impropers = Vec::with_capacity(ia.len() / 4);
    for i in 0..ia.len() / 4 {
        impropers.push(Improper {
            atoms: [
                at(ia, 4 * i)?,
                at(ia, 4 * i + 1)?,
                at(ia, 4 * i + 2)?,
                at(ia, 4 * i + 3)?,
            ],
            c0: at(ic, 4 * i)?,
            c1: at(ic, 4 * i + 1)?,
            c2: at(ic, 4 * i + 2)?,
            fc: at(ic, 4 * i + 3)?,
        });
    }

    // bonds (2 per) & angles (4 per: i, j=center, k, tripleFlag) for Stage-C constraints
    let bv = c.bonds();
    let mut bonds = Vec::with_capacity(bv.len() / 2);
    for i in 0..bv.len() / 2 {
        bonds.push([at(bv, 2 * i)?, at(bv, 2 * i + 1)?]);
    }
    let av = c.angles();
    let mut angles = Vec::with_capacity(av.len() / 4);
    for i in 0..av.len() / 4 {
        angles.push(Angle {
            atoms: [at(av, 4 * i)?, at(av, 4 * i + 1)?, at(av, 4 * i + 2)?],
            triple: at(av, 4 * i + 3)? != 0,
        });
    }

    // core-pin recipe: pinned atom set + conformer counts
    let pin_atoms: Vec<u32> = c.pin_atoms().iter().copied().collect();

    // double bonds for the acceptance checks
    let dbe = c.dbe();
    let mut double_bond_ends = Vec::with_capacity(dbe.len() / 3);
    for i in 0..dbe.len() / 3 {
        double_bond_ends.push([at(dbe, 3 * i)?, at(dbe, 3 * i + 1)?, at(dbe, 3 * i + 2)?]);
    }
    let (sa, ss) = (c.sdb_atoms(), c.sdb_sign());
    let mut stereo_double_bonds = Vec::with_capacity(ss.len());
    for i in 0..ss.len() {
        stereo_double_bonds.push(StereoDoubleBond {
            atoms: [
                at(sa, 4 * i)?,
                at(sa, 4 * i + 1)?,
                at(sa, 4 * i + 2)?,
                at(sa, 4 * i + 3)?,
            ],
            sign: at(ss, i)? as i8,
        });
    }

    Ok(MoleculeSpec {
        n_atoms: n,
        dim: 4,
        atomic_numbers,
        formal_charge: c.formal_charge(),
        bounds,
        // Native-only fields (f64 embed bounds; raw pre-smoothing matrix for the core-pin recipe).
        // This bridge spec is the oracle for gate comparison, never embedded, so they are left empty.
        bounds_f64: Vec::new(),
        raw_bounds: Vec::new(),
        raw_bounds_f64: Vec::new(),
        chiral_sets,
        tetrahedral_centers,
        exp_torsions,
        impropers,
        bonds,
        angles,
        bounds_force_scaling: 1.0, // EmbedParameters ETKDGv3 default
        pin_atoms,
        core_seeds: c.core_seeds(),
        sidechain_confs: c.sidechain_confs(),
        double_bond_ends,
        stereo_double_bonds,
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn heavy_atoms_via_rdkit_ffi() {
        assert_eq!(super::num_heavy_atoms("c1ccccc1"), 6);
        assert_eq!(super::num_heavy_atoms("not_a_smiles"), 0);
    }

    #[test]
    fn triangle_smooth_port_matches_rdkit() {
        // Validate bb-core's Rust triangle_smooth against RDKit's real triangleSmoothBounds on a
        // synthetic, triangle-satisfiable bounds matrix (loose bounds around fixed 3D points).
        let n = 8;
        let mut s = 0x1234_5678u64;
        let mut next = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 33) as f64 / (1u64 << 31) as f64 - 1.0) * 5.0
        };
        let pts: Vec<f64> = (0..n * 3).map(|_| next()).collect();
        let mut base = vec![0.0f64; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                let d = ((0..3)
                    .map(|c| (pts[i * 3 + c] - pts[j * 3 + c]).powi(2))
                    .sum::<f64>())
                .sqrt();
                base[i * n + j] = 1.5 * d; // upper
                base[j * n + i] = 0.5 * d; // lower
            }
        }
        for tol in [0.0, 0.05] {
            let rdkit = super::smooth_bounds(&base, n, tol).expect("RDKit smoothed the bounds");
            let mut mine: Vec<f32> = base.iter().map(|&x| x as f32).collect();
            assert!(bb_core::smooth::triangle_smooth(&mut mine, n, tol));
            let max_diff = (0..n * n)
                .map(|i| (mine[i] as f64 - rdkit[i]).abs())
                .fold(0.0f64, f64::max);
            assert!(
                max_diff < 1e-3,
                "tol {tol}: Rust port vs RDKit max diff {max_diff:.2e}"
            );
        }
    }

    #[test]
    fn spec_carries_torsions_and_chirality() {
        // cyclo(L-Ala)6 — 6 Cα stereocenters, ring amides (Divya's patched FC 80.0).
        let smi = "C[C@@H]1C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N1";
        let spec = super::build_spec(smi).expect("spec built");
        assert!(spec.n_atoms > 0);
        assert_eq!(spec.bounds.len(), spec.n_atoms * spec.n_atoms);
        assert!(spec.ub(0, 1) >= spec.lb(0, 1));

        // patched amide torsion FC 80.0 (proves Divya's engine via FFI, not stock 8.0)
        assert!(
            spec.exp_torsions
                .iter()
                .any(|t| t.v.iter().any(|&v| (v - 80.0).abs() < 0.1)),
            "patched amide torsion FC 80.0 not found"
        );

        // 6 tagged stereocenters (the 6 Cα); each chiral set is well-formed
        assert_eq!(
            spec.chiral_sets.len(),
            6,
            "expected 6 stereocenters for cyclo(L-Ala)6"
        );
        for cs in &spec.chiral_sets {
            assert!(cs.vol_hi > cs.vol_lo);
            assert!((cs.center as usize) < spec.n_atoms);
            assert!(cs.atoms.iter().all(|&a| (a as usize) < spec.n_atoms));
        }

        // UFF impropers on the sp2 amide C (=O) and N centers (3 permutations each)
        assert!(!spec.impropers.is_empty(), "no impropers");
        // carbonyl C (bound to sp2 O) -> fc 500/3 ≈ 166.7; N (not C-bound-to-O) -> fc 20
        assert!(
            spec.impropers
                .iter()
                .any(|im| (im.fc - 500.0 / 3.0).abs() < 1.0),
            "carbonyl-C improper fc (500/3) missing"
        );
        assert!(
            spec.impropers.iter().any(|im| (im.fc - 20.0).abs() < 0.5),
            "N improper fc (20) missing"
        );
    }
}
