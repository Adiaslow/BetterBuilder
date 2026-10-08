//! bb-core — the `MoleculeSpec` data contract: the distance-geometry problem for one molecule.
//!
//! Produced by setup — `bb-spec`, with no RDKit; `bb-rdkit`'s bridge builds the same contract for the
//! parity gates — and consumed by the `bb-embed` engine. Numeric arrays only; serde-serializable, so
//! the CLIs exchange it as JSON. [`smooth`] holds the triangle-inequality
//! operation on the bounds matrix — the one primitive both the spec builder and the engine share.

pub mod smooth;
pub mod vec3;

/// Signed chiral-volume constraint on four points. As an energy term (`vol_lo/vol_hi` non-zero) it
/// penalises the signed volume `(p_i−p_l)·((p_j−p_l)×(p_k−p_l))` outside `[vol_lo, vol_hi]`. With
/// `vol_lo == vol_hi == 0` it is an untagged tetrahedral center used only for the non-degeneracy
/// acceptance check.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ChiralSet {
    /// `d_idx0`: centroid, or the center atom itself for a 3-coordinate center.
    pub center: u32,
    /// The four substituent points defining the signed volume.
    pub atoms: [u32; 4],
    pub vol_lo: f64,
    pub vol_hi: f64,
    /// Scales the tetrahedral volume check ×0.25.
    pub fused_small_rings: bool,
}

/// CrystalFF M6 experimental-torsion term: `E = Σ_{n=0..5} V[n]·(1 + signs[n]·cos((n+1)·φ))`.
/// Atoms and coefficients come from RDKit `GetExperimentalTorsions`.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExpTorsion {
    pub atoms: [u32; 4],
    pub v: [f64; 6],
    pub signs: [i8; 6],
}

/// UFF improper / out-of-plane inversion term at an sp2 center (atom `j`, position 1). Energy
/// `E = fc·(C0 + C1·sinY + C2·cos2W)`. RDKit adds three per sp2-deg-3 C/N/O center (permutations);
/// `fc` already includes the DG force-scaling (×10). Coeffs from RDKit's UFF inversion function.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Improper {
    pub atoms: [u32; 4], // (i, j=center, k, l)
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub fc: f64,
}

/// A 1-3 angle (i-j-k, `j` central) from RDKit's `collectBondsAndAngles`. `triple` marks a
/// near-linear geometry (a triple bond, or two consecutive double bonds at a degree-2 center):
/// Stage C then constrains it to 179–180° rather than pinning the 1-3 distance.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Angle {
    pub atoms: [u32; 3], // (i, j=center, k)
    pub triple: bool,
}

/// A stereo double bond for the ETKDG acceptance check: controlling atoms `[stereoAtom0, begin,
/// end, stereoAtom1]` and `sign` (+1 = trans/E, −1 = cis/Z). The dihedral `stereoAtom0-begin-end-
/// stereoAtom1` must land on the correct side of 90° (`(dihedral − π/2)·sign ≥ 0`).
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct StereoDoubleBond {
    pub atoms: [u32; 4],
    pub sign: i8,
}

/// The per-molecule distance-geometry problem. Numeric arrays only. `bounds` is `n_atoms²`
/// row-major with the RDKit convention: **upper triangle (row<col) = upper bound, lower triangle
/// (row>col) = lower bound** (use [`MoleculeSpec::ub`] / [`MoleculeSpec::lb`]).
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct MoleculeSpec {
    pub n_atoms: usize,
    /// 3 or 4. Embedding is 4D when `useRandomCoords || n_chiral > 0`.
    pub dim: u8,
    /// Atomic number per atom, in the same order as every other index here (RDKit `SmilesToMol` +
    /// `AddHs`). Empty in specs written before this field existed.
    #[serde(default)]
    pub atomic_numbers: Vec<u8>,
    /// Net formal charge of the molecule.
    #[serde(default)]
    pub formal_charge: i32,
    pub bounds: Vec<f32>,
    /// The smoothed bounds matrix in f64 — the precision RDKit's force field actually uses. The embed
    /// reads these (via `ub64`/`lb64`) so its Stage-A/C distance and constraint terms match RDKit to
    /// machine precision, not the ~1e-5 an f32 bound amplifies to near a tight constraint. `bounds`
    /// stays f32 solely so the spec is byte-identical to the bridge for the spec gate. Empty in
    /// pre-field specs / the bridge oracle spec (never embedded).
    #[serde(default)]
    pub bounds_f64: Vec<f64>,
    /// The pre-smoothing `setTopolBounds` matrix (before triangle smoothing), `n*n` row-major.
    /// The core-pin recipe tightens *this* with the pinned distances and re-smooths at tol 0.05 —
    /// matching RDKit's `setupInitialBoundsMatrix` coordMap path, which also starts from the raw
    /// bounds (the tol-0.05 infeasibility repair is order-dependent, so starting from the already-
    /// smoothed `bounds` diverges). Empty in specs written before this field / for the free embed.
    #[serde(default)]
    pub raw_bounds: Vec<f32>,
    /// The pre-smoothing bounds matrix in f64 — what the core-pin recipe tightens (via
    /// `coord_map_bounds_f64`) so the sidechain bounds match RDKit's coordMap path to machine precision.
    /// `raw_bounds` stays f32 for symmetry with `bounds`; empty when unused.
    #[serde(default)]
    pub raw_bounds_f64: Vec<f64>,
    /// Tagged centers → chiral energy term.
    pub chiral_sets: Vec<ChiralSet>,
    /// Untagged C/N degree-4 centers → tetrahedral non-degeneracy check only.
    pub tetrahedral_centers: Vec<ChiralSet>,
    pub exp_torsions: Vec<ExpTorsion>,
    /// UFF impropers (sp2 planarity, Stage C).
    pub impropers: Vec<Improper>,
    /// Every bond (begin, end) over the AddHs mol. Stage C pins each 1-2 distance to the current
    /// (post-Stage-A) geometry ±0.01 Å; the long-range term uses these to know which pairs are 1-2.
    #[serde(default)]
    pub bonds: Vec<[u32; 2]>,
    /// Every 1-3 angle. Stage C pins each 1-3 distance to the current geometry (or 179–180° when
    /// `triple`); remaining pairs get bounds-matrix distance constraints.
    #[serde(default)]
    pub angles: Vec<Angle>,
    /// `EmbedParameters.boundsMatForceScaling` (ETKDGv3 = 1.0). Long-range Stage-C fc = ×10.
    #[serde(default = "default_bounds_force_scaling")]
    pub bounds_force_scaling: f32,
    /// Atoms pinned by the two-stage core-pin recipe: the largest-ring core + exocyclic `=O` and
    /// exocyclic C on ring N. Frozen to the seed conformer's coords while sidechains sample.
    #[serde(default)]
    pub pin_atoms: Vec<u32>,
    /// Number of core seed conformers (`rdkit_confs_1`): 20 if `exo≤2` else 10.
    #[serde(default)]
    pub core_seeds: u32,
    /// Sidechain conformers per core seed: 10 if `exo≤2` else 20. Total = `core_seeds × sidechain_confs`.
    #[serde(default)]
    pub sidechain_confs: u32,
    /// Double-bond substituent triples `(nbr, dbAtom, otherDbAtom)` — the linear-arrangement
    /// acceptance check (a substituent must not be collinear with the double bond).
    #[serde(default)]
    pub double_bond_ends: Vec<[u32; 3]>,
    /// Stereo double bonds — the cis/trans acceptance check.
    #[serde(default)]
    pub stereo_double_bonds: Vec<StereoDoubleBond>,
}

fn default_bounds_force_scaling() -> f32 {
    1.0
}

impl MoleculeSpec {
    /// Upper distance bound for pair (i, j). f32 storage, for the byte-identical spec gate.
    #[inline]
    pub fn ub(&self, i: usize, j: usize) -> f32 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.bounds[lo * self.n_atoms + hi]
    }

    /// Lower distance bound for pair (i, j). f32 storage, for the byte-identical spec gate.
    #[inline]
    pub fn lb(&self, i: usize, j: usize) -> f32 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.bounds[hi * self.n_atoms + lo]
    }

    /// Upper bound in the f64 precision RDKit's force field uses. Falls back to the f32 matrix when
    /// `bounds_f64` is absent (older specs / the bridge oracle spec, which is never embedded).
    #[inline]
    pub fn ub64(&self, i: usize, j: usize) -> f64 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        if self.bounds_f64.is_empty() {
            self.bounds[lo * self.n_atoms + hi] as f64
        } else {
            self.bounds_f64[lo * self.n_atoms + hi]
        }
    }

    /// Lower bound in f64 precision (see [`Self::ub64`]).
    #[inline]
    pub fn lb64(&self, i: usize, j: usize) -> f64 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        if self.bounds_f64.is_empty() {
            self.bounds[hi * self.n_atoms + lo] as f64
        } else {
            self.bounds_f64[hi * self.n_atoms + lo]
        }
    }

    /// A per-seed variant for the core-pin recipe: identical topology, distance bounds replaced by the
    /// coordMap-tightened `bounds_f64` (with its f32 shadow, which the acceptance checks read via
    /// [`Self::ub`]/[`Self::lb`]). The pre-smoothing `raw_bounds`/`raw_bounds_f64` are left empty — the
    /// sidechain embed never reads them — so this avoids the two n² matrices `Clone` would copy per
    /// seed only to discard. Full struct literal on purpose: a new field forces a decision here.
    pub fn with_seed_bounds(&self, bounds_f64: Vec<f64>) -> MoleculeSpec {
        MoleculeSpec {
            bounds: bounds_f64.iter().map(|&x| x as f32).collect(),
            bounds_f64,
            raw_bounds: Vec::new(),
            raw_bounds_f64: Vec::new(),
            n_atoms: self.n_atoms,
            dim: self.dim,
            atomic_numbers: self.atomic_numbers.clone(),
            formal_charge: self.formal_charge,
            chiral_sets: self.chiral_sets.clone(),
            tetrahedral_centers: self.tetrahedral_centers.clone(),
            exp_torsions: self.exp_torsions.clone(),
            impropers: self.impropers.clone(),
            bonds: self.bonds.clone(),
            angles: self.angles.clone(),
            bounds_force_scaling: self.bounds_force_scaling,
            pin_atoms: self.pin_atoms.clone(),
            core_seeds: self.core_seeds,
            sidechain_confs: self.sidechain_confs,
            double_bond_ends: self.double_bond_ends.clone(),
            stereo_double_bonds: self.stereo_double_bonds.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_triangle_convention() {
        // 2 atoms; UB(0,1) in upper-tri [0*2+1], LB(0,1) in lower-tri [1*2+0].
        let spec = MoleculeSpec {
            n_atoms: 2,
            dim: 4,
            bounds: vec![0.0, 1.46 /*UB*/, 1.44 /*LB*/, 0.0],
            ..Default::default()
        };
        assert_eq!(spec.ub(0, 1), 1.46);
        assert_eq!(spec.lb(0, 1), 1.44);
        assert_eq!(spec.ub(1, 0), 1.46); // order-independent
        assert!(spec.ub(0, 1) >= spec.lb(0, 1));
    }

    /// A spec survives a JSON round trip bit for bit, its f64 matrices included: the `MoleculeSpec`
    /// JSON contract between `bb-spec-native` and `bb-embed` must hand the embed exactly the spec that
    /// was built, or the chaotic embed diverges from the in-process pipeline. The values are full-
    /// precision f64s spread over the range bounds take, 1 to 20 Å.
    #[test]
    fn molecule_spec_round_trips_json() {
        let n = 100;
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            1.0 + 19.0 * ((state >> 11) as f64 / (1u64 << 53) as f64)
        };
        let bounds_f64: Vec<f64> = (0..n * n).map(|_| next()).collect();
        let raw_bounds_f64: Vec<f64> = (0..n * n).map(|_| next()).collect();
        let spec = MoleculeSpec {
            n_atoms: n,
            dim: 4,
            bounds: bounds_f64.iter().map(|&x| x as f32).collect(),
            bounds_f64,
            raw_bounds_f64,
            ..Default::default()
        };
        let back: MoleculeSpec = serde_json::from_str(&serde_json::to_string(&spec).unwrap()).unwrap();
        let changed = |a: &[f64], b: &[f64]| a.iter().zip(b).filter(|(x, y)| x.to_bits() != y.to_bits()).count();
        assert_eq!(back.n_atoms, n);
        assert_eq!(changed(&back.bounds_f64, &spec.bounds_f64), 0, "bounds_f64 values changed");
        assert_eq!(changed(&back.raw_bounds_f64, &spec.raw_bounds_f64), 0, "raw_bounds_f64 values changed");
        assert_eq!(back.bounds, spec.bounds, "f32 bounds changed");
    }
}
