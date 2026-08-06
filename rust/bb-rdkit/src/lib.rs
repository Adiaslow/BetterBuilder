//! bb-rdkit — FFI to the patched RDKit C++ (2026.09.1pre + Divya's amide torsions) for setup.
//!
//! The candidate is pure Rust with zero Python; setup/perception is best served by RDKit, so we
//! bind its C++ directly with `cxx` (best-in-class Rust↔C++ interop). This turns a SMILES into a
//! full [`bb_core::MoleculeSpec`] — bounds matrix + experimental torsions (which carry Divya's
//! amide-planarity reweighting) + chiral/tetrahedral sets + UFF impropers + bonds/angles +
//! double-bond stereo + the core-pin recipe counts — that the `bb-embed` minimizer consumes. Links
//! our *from-source patched* RDKit build.

use bb_core::{Angle, ChiralSet, ExpTorsion, Improper, MoleculeSpec, StereoDoubleBond};
use cxx::vector::VectorElement;
use cxx::CxxVector;
use thiserror::Error;

/// Errors from setup extraction. The only user-recoverable case is [`SpecError::InvalidSmiles`];
/// [`SpecError::MalformedFfi`] is defensive — it fires only if the C++ bridge returns an array whose
/// length is not the expected multiple, which would be a bug in `bridge.cc`, surfaced (not panicked)
/// so a library caller (e.g. the Python binding) gets an error instead of a crashed process.
#[derive(Debug, Error)]
pub enum SpecError {
    #[error("RDKit could not parse SMILES: {0:?}")]
    InvalidSmiles(String),
    #[error("malformed data from the RDKit FFI bridge: index {0} out of range")]
    MalformedFfi(usize),
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

        fn num_heavy_atoms(smiles: &str) -> usize;

        /// Test-only: RDKit's real `triangleSmoothBounds` on a flat n×n bounds matrix.
        fn smooth_bounds(bounds: &[f64], n: usize, tol: f64) -> Vec<f64>;

        /// Opaque holder for one molecule's extracted distance-geometry problem.
        type MolSpecC;
        fn build_mol_spec(smiles: &str) -> UniquePtr<MolSpecC>;
        fn ok(self: &MolSpecC) -> bool;
        fn n_atoms(self: &MolSpecC) -> usize;
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

/// Heavy-atom count via RDKit's SMILES parser (the FFI link check). 0 on parse failure.
pub fn num_heavy_atoms(smiles: &str) -> usize {
    ffi::num_heavy_atoms(smiles)
}

/// RDKit's real `triangleSmoothBounds` on a flat `n×n` bounds matrix (row<col = upper, row>col =
/// lower). Used to validate bb-embed's Rust port; not on the embed path.
pub fn smooth_bounds(bounds: &[f64], n: usize, tol: f64) -> Vec<f64> {
    ffi::smooth_bounds(bounds, n, tol)
}

/// Build a [`MoleculeSpec`] from a SMILES using the patched RDKit: bounds matrix + experimental
/// torsions + chiral/tetrahedral sets + UFF impropers + bonds/angles + double-bond stereo. Returns
/// [`SpecError::InvalidSmiles`] if RDKit cannot parse the SMILES.
pub fn build_spec(smiles: &str) -> Result<MoleculeSpec, SpecError> {
    let c = ffi::build_mol_spec(smiles);
    if !c.ok() {
        return Err(SpecError::InvalidSmiles(smiles.to_string()));
    }
    let n = c.n_atoms();
    let bounds: Vec<f32> = c.bounds().iter().map(|&x| x as f32).collect();

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
        let mut v = [0f32; 6];
        let mut signs = [0i8; 6];
        for k in 0..6 {
            v[k] = at(tv, 6 * i + k)? as f32;
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
            vol_lo: at(cv, 2 * i)? as f32,
            vol_hi: at(cv, 2 * i + 1)? as f32,
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
            c0: at(ic, 4 * i)? as f32,
            c1: at(ic, 4 * i + 1)? as f32,
            c2: at(ic, 4 * i + 2)? as f32,
            fc: at(ic, 4 * i + 3)? as f32,
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
        bounds,
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
        // Validate bb-embed's Rust triangle_smooth against RDKit's real triangleSmoothBounds on a
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
            let rdkit = super::smooth_bounds(&base, n, tol);
            let mut mine: Vec<f32> = base.iter().map(|&x| x as f32).collect();
            assert!(bb_embed::bounds::triangle_smooth(&mut mine, n, tol));
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
