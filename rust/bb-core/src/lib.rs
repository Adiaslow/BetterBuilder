//! bb-core — the `MoleculeSpec` contract: the distance-geometry problem for one molecule.
//!
//! Produced from RDKit's perception (bounds matrix + experimental torsions + chirality) and
//! consumed by the `bb-embed` CPU ETKDG minimizer (the BetterBuilder candidate). RDKit does all
//! perception; this side is pure numeric geometry. Serde-serializable so setup (Python/RDKit) and
//! the embed (Rust) exchange it as JSON.

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
    pub vol_lo: f32,
    pub vol_hi: f32,
    /// Scales the tetrahedral volume check ×0.25.
    pub fused_small_rings: bool,
}

/// CrystalFF M6 experimental-torsion term: `E = Σ_{n=0..5} V[n]·(1 + signs[n]·cos((n+1)·φ))`.
/// Atoms + coefficients come straight from RDKit `GetExperimentalTorsions`.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExpTorsion {
    pub atoms: [u32; 4],
    pub v: [f32; 6],
    pub signs: [i8; 6],
}

/// UFF improper / out-of-plane inversion term at an sp2 center (atom `j`, position 1). Energy
/// `E = fc·(C0 + C1·sinY + C2·cos2W)`. RDKit adds three per sp2-deg-3 C/N/O center (permutations);
/// `fc` already includes the DG force-scaling (×10). Coeffs from RDKit's UFF inversion function.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Improper {
    pub atoms: [u32; 4], // (i, j=center, k, l)
    pub c0: f32,
    pub c1: f32,
    pub c2: f32,
    pub fc: f32,
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
    pub bounds: Vec<f32>,
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
    /// Upper distance bound for pair (i, j).
    #[inline]
    pub fn ub(&self, i: usize, j: usize) -> f32 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.bounds[lo * self.n_atoms + hi]
    }

    /// Lower distance bound for pair (i, j).
    #[inline]
    pub fn lb(&self, i: usize, j: usize) -> f32 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.bounds[hi * self.n_atoms + lo]
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

    #[test]
    fn molecule_spec_round_trips_json() {
        let spec = MoleculeSpec {
            n_atoms: 2,
            dim: 4,
            bounds: vec![0.0, 1.46, 1.44, 0.0],
            ..Default::default()
        };
        let json = serde_json::to_string(&spec).unwrap();
        let back: MoleculeSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n_atoms, 2);
        assert_eq!(back.ub(0, 1), 1.46);
    }
}
