//! Force-field terms for the staged embed, in f64 with analytic gradients.
//!
//! Stages A and B ([`stage_a_energy_grad`], [`stage_b_energy_grad`]): distance-bound violation over
//! all dims + chiral volume (x,y,z) + 4th-dimension penalty, differing only in the chiral and
//! 4th-dimension weights. Stage C ([`stage_c_energy_grad`]): M6 experimental torsions + UFF
//! impropers + flat-bottom distance and angle constraints.
//!
//! The chiral and 4th-dimension terms use energy `0.5·w·x²` with gradient `w·x`, as in RDKit's
//! `ChiralViolationContribs` / `FourthDimContribs`.

use std::sync::atomic::{AtomicU64, Ordering};

use bb_core::vec3::{cross, dot as dot3, norm as norm3, scale as scale3, sub as subv};
use bb_core::MoleculeSpec;

/// Count of DistGeom (Stage A/B) and Stage-C force-field evaluations. Only incremented when the
/// `profile` feature is enabled.
pub static FF_EVALS: AtomicU64 = AtomicU64::new(0);

/// Count one force-field evaluation (no-op unless the `profile` feature is on).
#[inline]
fn ff_tick() {
    #[cfg(feature = "profile")]
    FF_EVALS.fetch_add(1, Ordering::Relaxed);
}
#[inline]
pub fn ff_evals() -> u64 {
    FF_EVALS.load(Ordering::Relaxed)
}
#[inline]
pub fn reset_ff_evals() {
    FF_EVALS.store(0, Ordering::Relaxed);
}

pub const W_DIST: f64 = 1.0;
pub const W_CHIRAL: f64 = 1.0;
pub const W_FOURTH: f64 = 0.1;

/// Difference of two atoms read from a strided coordinate slice: `p[a] − p[b]` (x,y,z only).
#[inline]
fn sub3(p: &[f64], a: usize, b: usize, dim: usize) -> [f64; 3] {
    [
        p[a * dim] - p[b * dim],
        p[a * dim + 1] - p[b * dim + 1],
        p[a * dim + 2] - p[b * dim + 2],
    ]
}

/// Basin threshold (`basinSizeTol`): distance terms are dropped when `ub−lb` exceeds it. RDKit uses
/// 5.0 by default and 1e8 (i.e. keep everything) when `useRandomCoords`.
pub const BASIN_DEFAULT: f64 = 5.0;
pub const BASIN_ALL: f64 = 1e8;

/// Stage-A DistGeom objective: distance-bound violation + chiral (weight 1.0) + 4th-dim (weight 0.1).
/// (`constructForceField(mmat, pos, csets, 1.0, 0.1, …, basinSizeTol)`.)
pub fn stage_a_energy_grad(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    basin: f64,
) -> (f64, Vec<f64>) {
    dist_geom_energy_grad(spec, coords, dim, W_CHIRAL, W_FOURTH, basin)
}

/// Native's chiral + 4th-dimension Stage-A energy in native's own (gradient-consistent, `½·w·x²`)
/// convention, at Stage-A weights. This is exactly the amount by which RDKit's `calcEnergy` exceeds
/// native's total energy (RDKit's Chiral/FourthDim `getEnergy` use `w·x²`, i.e. twice native's) — see
/// [`stage_a_reject_energy`].
fn chiral_fourth_energy(spec: &MoleculeSpec, coords: &[f64], dim: usize) -> f64 {
    let mut e = 0.0;
    for cs in &spec.chiral_sets {
        let (vol, bound, _) = chiral_vol(coords, cs, dim);
        if let Some(b) = bound {
            e += 0.5 * W_CHIRAL * (vol - b) * (vol - b);
        }
    }
    if dim == 4 {
        for i in 0..spec.n_atoms {
            let w = coords[i * dim + 3];
            e += 0.5 * W_FOURTH * w * w;
        }
    }
    e
}

/// The signed chiral volume of one center, its clamped bound (`Some(limit)` when `vol` is outside
/// `[vol_lo, vol_hi]`, else `None`), and the three edge vectors `[i−l, j−l, k−l]` the gradient reuses.
/// Single source for the chiral-volume test shared by the gradient accumulation and the energy-only
/// reject path.
#[inline]
fn chiral_vol(coords: &[f64], cs: &bb_core::ChiralSet, dim: usize) -> (f64, Option<f64>, [[f64; 3]; 3]) {
    let (i, j, k, l) = (
        cs.atoms[0] as usize,
        cs.atoms[1] as usize,
        cs.atoms[2] as usize,
        cs.atoms[3] as usize,
    );
    let v1 = sub3(coords, i, l, dim);
    let v2 = sub3(coords, j, l, dim);
    let v3 = sub3(coords, k, l, dim);
    let vol = dot3(v1, cross(v2, v3));
    let bound = if vol < cs.vol_lo {
        Some(cs.vol_lo)
    } else if vol > cs.vol_hi {
        Some(cs.vol_hi)
    } else {
        None
    };
    (vol, bound, [v1, v2, v3])
}

/// RDKit's Stage-A `calcEnergy` — the value `firstMinimization` thresholds against
/// `MAX_MINIMIZED_E_PER_ATOM` (0.05/atom) for its per-atom energy reject. RDKit's Chiral/FourthDim
/// `getEnergy` use `w·x²` (no ½), *inconsistent* with their `w·x` gradient; native's own energy uses
/// the gradient-consistent `½·w·x²`. So RDKit's reject energy is native's total energy plus the
/// chiral+4th contribution once more (the distance term matches). Used ONLY to reproduce RDKit's
/// reject decision faithfully — never as an optimization objective (native minimizes its own,
/// gradient-consistent energy). Verified boundary-exact against RDKit's `calcEnergy` bridge.
pub fn stage_a_reject_energy(spec: &MoleculeSpec, coords: &[f64], dim: usize, basin: f64) -> f64 {
    stage_a_energy_grad(spec, coords, dim, basin).0 + chiral_fourth_energy(spec, coords, dim)
}

/// Distance-bound violation term over all dims, accumulated into `e`/`g` for every atom pair. Pairs
/// whose bounds are wider than `basin` are skipped.
#[inline]
fn accum_dist_term(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    basin: f64,
    e: &mut f64,
    g: &mut [f64],
) {
    let n = spec.n_atoms;
    for i in 0..n {
        for j in (i + 1)..n {
            let ub = spec.ub64(i, j);
            let lb = spec.lb64(i, j);
            if ub - lb > basin {
                continue; // basin filter: skip loosely-bounded (flexible) pairs
            }
            accum_one_pair(coords, dim, i, j, lb * lb, ub * ub, e, g);
        }
    }
}

/// A distance-bound atom pair that survives the basin filter: `(i, j, lb², ub²)`. The surviving set
/// and its squared bounds are constant across a whole minimization, so [`build_dist_pairs`] computes
/// it once and the hot Stage-A/B eval iterates it instead of re-scanning all n²/2 pairs per call.
pub type DistPair = (usize, usize, f64, f64);

/// Basin-surviving distance pairs for `basin`, in the exact `(i<j)` row-major order the per-eval scan
/// visited — so accumulating over them reproduces the scan's summation order bit-for-bit.
pub fn build_dist_pairs(spec: &MoleculeSpec, basin: f64) -> Vec<DistPair> {
    let n = spec.n_atoms;
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let ub = spec.ub64(i, j);
            let lb = spec.lb64(i, j);
            if ub - lb > basin {
                continue;
            }
            pairs.push((i, j, lb * lb, ub * ub));
        }
    }
    pairs
}

/// Distance-bound contribution of one pair with squared bounds `(lb2, ub2)`, accumulated into `e`/`g`.
/// The single source for the Stage-A/B distance term: the basin scan and the pre-built pair list both
/// dispatch here, so the physics — and its exact FP rounding — exist in one place. `d = √d2` is taken
/// once per violating pair and reused for both `pre` and the gradient normalization.
#[inline]
fn accum_one_pair(
    coords: &[f64],
    dim: usize,
    i: usize,
    j: usize,
    lb2: f64,
    ub2: f64,
    e: &mut f64,
    g: &mut [f64],
) {
    let mut d2 = 0.0;
    for c in 0..dim {
        let d = coords[i * dim + c] - coords[j * dim + c];
        d2 += d * d;
    }
    let (val, pre, d) = if d2 > ub2 {
        let d = d2.sqrt();
        (d2 / ub2 - 1.0, 4.0 * (d2 / ub2 - 1.0) * (d / ub2), d)
    } else if d2 < lb2 {
        let d = d2.sqrt();
        let s = d2 + lb2;
        (
            2.0 * lb2 / s - 1.0,
            8.0 * lb2 * d * (1.0 - 2.0 * lb2 / s) / (s * s),
            d,
        )
    } else {
        return; // within bounds: no contribution (matches the old (0,0) + `val > 0` guard)
    };
    if val > 0.0 {
        *e += W_DIST * val * val;
        if d > 1e-12 {
            for c in 0..dim {
                let gg = W_DIST * pre * (coords[i * dim + c] - coords[j * dim + c]) / d;
                g[i * dim + c] += gg;
                g[j * dim + c] -= gg;
            }
        }
    }
}

/// Stage-A/B distance term over a pre-built basin-surviving pair list (the hot-path form of
/// [`accum_dist_term`]).
#[inline]
fn accum_dist_term_pairs(coords: &[f64], dim: usize, pairs: &[DistPair], e: &mut f64, g: &mut [f64]) {
    for &(i, j, lb2, ub2) in pairs {
        accum_one_pair(coords, dim, i, j, lb2, ub2, e, g);
    }
}

/// Chiral-volume (weight `w_chiral`, x/y/z) and fourth-dimension (weight `w_fourth`, 4th coord)
/// penalties — the part of the DistGeom objective shared by the basin-scan and pre-built-pair evals.
#[inline]
fn accum_chiral_fourth(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    w_chiral: f64,
    w_fourth: f64,
    e: &mut f64,
    g: &mut [f64],
) {
    // --- chiral-volume violation (x,y,z only) ---
    for cs in &spec.chiral_sets {
        let (vol, bound, [v1, v2, v3]) = chiral_vol(coords, cs, dim);
        if let Some(b) = bound {
            let (i, j, k, l) = (
                cs.atoms[0] as usize,
                cs.atoms[1] as usize,
                cs.atoms[2] as usize,
                cs.atoms[3] as usize,
            );
            *e += 0.5 * w_chiral * (vol - b) * (vol - b);
            let pre = w_chiral * (vol - b); // matches RDKit ChiralViolationContribs (no factor of 2)
            let d_i = cross(v2, v3);
            let d_j = cross(v3, v1);
            let d_k = cross(v1, v2);
            for c in 0..3 {
                g[i * dim + c] += pre * d_i[c];
                g[j * dim + c] += pre * d_j[c];
                g[k * dim + c] += pre * d_k[c];
                g[l * dim + c] += pre * (-d_i[c] - d_j[c] - d_k[c]);
            }
        }
    }

    // --- fourth-dimension penalty (4th coord only) ---
    if dim == 4 {
        for i in 0..spec.n_atoms {
            let w = coords[i * dim + 3];
            *e += 0.5 * w_fourth * w * w;
            g[i * dim + 3] += w_fourth * w; // matches RDKit FourthDimContribs (no factor of 2)
        }
    }
}

/// Stage-B DistGeom objective (`minimizeFourthDimension`): the same terms with the chiral weight
/// dropped to 0.2 and the 4th-dim weight raised to 1.0, squeezing out the 4th dimension.
/// (`constructForceField(mmat, pos, csets, 0.2, 1.0, …, basinSizeTol)`.)
pub fn stage_b_energy_grad(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    basin: f64,
) -> (f64, Vec<f64>) {
    dist_geom_energy_grad(spec, coords, dim, 0.2, 1.0, basin)
}

/// DistGeom energy and gradient with tunable chiral / 4th-dim weights and basin threshold. `coords`
/// is `n_atoms * dim` row-major, `dim` ∈ {3,4}.
fn dist_geom_energy_grad(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    w_chiral: f64,
    w_fourth: f64,
    basin: f64,
) -> (f64, Vec<f64>) {
    ff_tick();
    let n = spec.n_atoms;
    debug_assert_eq!(coords.len(), n * dim);
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];

    // --- distance-bound violation (all dims), weight 1.0 ---
    accum_dist_term(spec, coords, dim, basin, &mut e, &mut g);
    // --- chiral-volume (x,y,z) + fourth-dimension (4th coord) penalties ---
    accum_chiral_fourth(spec, coords, dim, w_chiral, w_fourth, &mut e, &mut g);

    (e, g)
}

/// Same DistGeom objective as [`dist_geom_energy_grad`] but iterating a pre-built basin-surviving
/// [`DistPair`] list instead of re-scanning all pairs — the hot-path form used by the minimizer, which
/// builds the list once per Stage-A/B minimization. Bit-identical to the scan form (same surviving
/// pairs, same order, same physics via [`accum_one_pair`]).
fn dist_geom_energy_grad_pairs(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    w_chiral: f64,
    w_fourth: f64,
    pairs: &[DistPair],
) -> (f64, Vec<f64>) {
    ff_tick();
    debug_assert_eq!(coords.len(), spec.n_atoms * dim);
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];
    accum_dist_term_pairs(coords, dim, pairs, &mut e, &mut g);
    accum_chiral_fourth(spec, coords, dim, w_chiral, w_fourth, &mut e, &mut g);
    (e, g)
}

/// Stage-A objective over a pre-built [`DistPair`] list (see [`stage_a_energy_grad`]).
pub fn stage_a_energy_grad_pairs(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    pairs: &[DistPair],
) -> (f64, Vec<f64>) {
    dist_geom_energy_grad_pairs(spec, coords, dim, W_CHIRAL, W_FOURTH, pairs)
}

/// Stage-B objective over a pre-built [`DistPair`] list (see [`stage_b_energy_grad`]).
pub fn stage_b_energy_grad_pairs(
    spec: &MoleculeSpec,
    coords: &[f64],
    dim: usize,
    pairs: &[DistPair],
) -> (f64, Vec<f64>) {
    dist_geom_energy_grad_pairs(spec, coords, dim, 0.2, 1.0, pairs)
}

/// A point read from a strided coordinate slice (x,y,z only).
#[inline]
fn pt(p: &[f64], i: usize, dim: usize) -> [f64; 3] {
    [p[i * dim], p[i * dim + 1], p[i * dim + 2]]
}

/// M6 experimental-torsion term over `spec.exp_torsions` (x,y,z only):
/// `E = Σ_{m=1..6} V[m-1]·(1 + signs[m-1]·cos(m·φ))`, φ = signed dihedral(i,j,k,l). Torsions whose
/// geometry is collinear are skipped.
pub fn torsion_energy_grad(spec: &MoleculeSpec, coords: &[f64], dim: usize) -> (f64, Vec<f64>) {
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];
    // RDKit `TorsionAngleContribs::getGrad` *returns* from the whole loop on the first torsion whose
    // plane-normal magnitude is `isDoubleZero` (collinear, e.g. through an sp nitrile/alkyne),
    // abandoning the gradient of that torsion AND every torsion after it in emission order. Its
    // *energy* (`calcTorsionEnergyM6`) has no such early-out. We reproduce that behavior exactly:
    // `grad_dead` latches at the first degenerate torsion and suppresses only the gradient thereafter,
    // while the energy keeps accumulating. This depends on `spec.exp_torsions` being in RDKit's
    // emission order — which native's assignment reproduces by construction (library file order ×
    // match order × first-claim dedup); the corpus parity gate is the check on that.
    let mut grad_dead = false;
    for t in &spec.exp_torsions {
        let (i, j, k, l) = (
            t.atoms[0] as usize,
            t.atoms[1] as usize,
            t.atoms[2] as usize,
            t.atoms[3] as usize,
        );
        let (pi, pj, pk, pl) = (
            pt(coords, i, dim),
            pt(coords, j, dim),
            pt(coords, k, dim),
            pt(coords, l, dim),
        );
        // RDKit `TorsionAngleContribs::getGrad` convention: r0=i-j, r1=k-j, r2=j-k, r3=l-k;
        // the two plane normals are t0=r0×r1, t1=r2×r3 (equal to native's ∓n1,∓n2, so cosφ agrees).
        let r0 = subv(pi, pj);
        let r1 = subv(pk, pj);
        let r2 = subv(pj, pk);
        let r3 = subv(pl, pk);
        let t0 = cross(r0, r1);
        let t1 = cross(r2, r3);
        let (d0, d1) = (norm3(t0), norm3(t1));
        if d0 < 1e-10 || d1 < 1e-10 {
            // RDKit returns here: this torsion and all after it contribute no gradient.
            grad_dead = true;
            continue; // degenerate (collinear) — no energy or gradient for this one either
        }
        let t0n = scale3(t0, 1.0 / d0);
        let t1n = scale3(t1, 1.0 / d1);
        let cos_phi = dot3(t0n, t1n).clamp(-1.0, 1.0);
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        let sin_phi = if sin_phi_sq > 0.0 { sin_phi_sq.sqrt() } else { 0.0 };
        let phi = cos_phi.acos();

        // Energy: E = Σ V[m-1]·(1 + signs[m-1]·cos(m·φ)). cos is even, so |φ| suffices.
        for idx in 0..6 {
            let m = (idx + 1) as f64;
            e += t.v[idx] * (1.0 + t.signs[idx] as f64 * (m * phi).cos());
        }

        // gradient abandoned once RDKit's getGrad has returned (at an earlier degenerate torsion)
        if grad_dead {
            continue;
        }

        // dE/dφ via RDKit's exact Chebyshev-derivative polynomial. NB: RDKit's m=6 term reuses
        // forceConstants[4]/signs[4] (not [5]) — an energy/gradient inconsistency in RDKit itself;
        // we replicate it so the gradient (what the minimizer follows) matches bit-for-bit.
        let (c1, c2, c3, c4, c5) = (
            cos_phi,
            cos_phi.powi(2),
            cos_phi.powi(3),
            cos_phi.powi(4),
            cos_phi.powi(5),
        );
        let (fc, sg) = (|x: usize| t.v[x], |x: usize| t.signs[x] as f64);
        let de_dphi = -fc(0) * sg(0) * sin_phi
            - 2.0 * fc(1) * sg(1) * (2.0 * c1 * sin_phi)
            - 3.0 * fc(2) * sg(2) * (4.0 * c2 * sin_phi - sin_phi)
            - 4.0 * fc(3) * sg(3) * (8.0 * c3 * sin_phi - 4.0 * c1 * sin_phi)
            - 5.0 * fc(4) * sg(4) * (16.0 * c4 * sin_phi - 12.0 * c2 * sin_phi + sin_phi)
            - 6.0 * fc(4) * sg(4) * (32.0 * c5 * sin_phi - 32.0 * c3 * sin_phi + 6.0 * sin_phi);

        // sinTerm = -dE/dφ · 1/sinφ analytically cancels the 1/sinφ in dφ/dx (the source of the old
        // form's instability near collinear); near sinφ=0 RDKit substitutes 1/cosφ (Niketic & Rasmussen).
        let sin_term = -de_dphi * if sin_phi < 1e-10 { 1.0 / cos_phi } else { 1.0 / sin_phi };

        // calcTorsionGrad (UFF/TorsionAngle.cpp): d(cosφ)/dt projected onto the four atoms.
        let dct = [
            (t1n[0] - cos_phi * t0n[0]) / d0,
            (t1n[1] - cos_phi * t0n[1]) / d0,
            (t1n[2] - cos_phi * t0n[2]) / d0,
            (t0n[0] - cos_phi * t1n[0]) / d1,
            (t0n[1] - cos_phi * t1n[1]) / d1,
            (t0n[2] - cos_phi * t1n[2]) / d1,
        ];
        let gi = [
            sin_term * (dct[2] * r1[1] - dct[1] * r1[2]),
            sin_term * (dct[0] * r1[2] - dct[2] * r1[0]),
            sin_term * (dct[1] * r1[0] - dct[0] * r1[1]),
        ];
        let gj = [
            sin_term * (dct[1] * (r1[2] - r0[2]) + dct[2] * (r0[1] - r1[1]) - dct[4] * r3[2] + dct[5] * r3[1]),
            sin_term * (dct[0] * (r0[2] - r1[2]) + dct[2] * (r1[0] - r0[0]) + dct[3] * r3[2] - dct[5] * r3[0]),
            sin_term * (dct[0] * (r1[1] - r0[1]) + dct[1] * (r0[0] - r1[0]) - dct[3] * r3[1] + dct[4] * r3[0]),
        ];
        let gk = [
            sin_term * (dct[1] * r0[2] - dct[2] * r0[1] + dct[4] * (r3[2] - r2[2]) + dct[5] * (r2[1] - r3[1])),
            sin_term * (-dct[0] * r0[2] + dct[2] * r0[0] + dct[3] * (r2[2] - r3[2]) + dct[5] * (r3[0] - r2[0])),
            sin_term * (dct[0] * r0[1] - dct[1] * r0[0] + dct[3] * (r3[1] - r2[1]) + dct[4] * (r2[0] - r3[0])),
        ];
        let gl = [
            sin_term * (dct[4] * r2[2] - dct[5] * r2[1]),
            sin_term * (dct[5] * r2[0] - dct[3] * r2[2]),
            sin_term * (dct[3] * r2[1] - dct[4] * r2[0]),
        ];
        for c in 0..3 {
            g[i * dim + c] += gi[c];
            g[j * dim + c] += gj[c];
            g[k * dim + c] += gk[c];
            g[l * dim + c] += gl[c];
        }
    }
    (e, g)
}

/// UFF out-of-plane / inversion term, the Stage-C planarity contribution (RDKit
/// `InversionContribs`). Per contrib on atoms (i, j=center, k, l):
/// `E = fc·(C0 + C1·sinY + C2·cos2W)`, where `sinY` is the sine of the Wilson out-of-plane angle
/// (`cosY = n̂·r̂JL`, `n̂ ⟂` the i-j-k plane) and `cos2W = 2·sinY² − 1`. Coefficients come from
/// `spec.impropers`.
///
/// Near-degenerate geometries (sinY→0 or sinθ→0) are handled as RDKit does — sinY and sinθ are
/// floored at 1e-8 in the gradient and the large-but-finite value is kept, not skipped — so the
/// force the minimizer follows matches RDKit's bit-for-bit even at singular conformers.
/// Per-improper geometry shared by the energy and gradient: unit bond vectors from the center `j`,
/// their lengths, and the out-of-plane angle `cosY`/`sin²Y`. `None` for a degenerate geometry (a
/// zero-length bond or a collinear i-j-k giving a zero plane normal) that contributes nothing — the
/// caller skips it, matching the original `continue`s.
struct ImproperGeom {
    rji: [f64; 3],
    rjk: [f64; 3],
    rjl: [f64; 3],
    dji: f64,
    djk: f64,
    djl: f64,
    cosy: f64,
    siny_sq: f64,
}

#[inline]
fn improper_geom(coords: &[f64], imp: &bb_core::Improper, dim: usize) -> Option<ImproperGeom> {
    let (i, j, k, l) = (
        imp.atoms[0] as usize,
        imp.atoms[1] as usize,
        imp.atoms[2] as usize,
        imp.atoms[3] as usize,
    );
    let (p1, p2, p3, p4) = (
        pt(coords, i, dim),
        pt(coords, j, dim),
        pt(coords, k, dim),
        pt(coords, l, dim),
    );
    // bond vectors from the center j
    let mut rji = subv(p1, p2);
    let mut rjk = subv(p3, p2);
    let mut rjl = subv(p4, p2);
    let (dji, djk, djl) = (norm3(rji), norm3(rjk), norm3(rjl));
    if dji < 1e-8 || djk < 1e-8 || djl < 1e-8 {
        return None;
    }
    rji = scale3(rji, 1.0 / dji);
    rjk = scale3(rjk, 1.0 / djk);
    rjl = scale3(rjl, 1.0 / djl);
    // plane normal n̂ = (-r̂JI) × r̂JK
    let mut n = cross(scale3(rji, -1.0), rjk);
    let nn = norm3(n);
    if nn < 1e-8 {
        return None;
    }
    n = scale3(n, 1.0 / nn);
    let cosy = dot3(n, rjl).clamp(-1.0, 1.0);
    let siny_sq = 1.0 - cosy * cosy;
    Some(ImproperGeom { rji, rjk, rjl, dji, djk, djl, cosy, siny_sq })
}

/// Energy of one improper from its [`ImproperGeom`] (RDKit: `sinY = sqrt(sin²Y)` if `>0` else `0`).
#[inline]
fn improper_energy_one(imp: &bb_core::Improper, gm: &ImproperGeom) -> f64 {
    let siny_e = if gm.siny_sq > 0.0 { gm.siny_sq.sqrt() } else { 0.0 };
    let cos2w = 2.0 * siny_e * siny_e - 1.0;
    imp.fc * (imp.c0 + imp.c1 * siny_e + imp.c2 * cos2w)
}

/// UFF improper (sp2 planarity) energy only — no gradient allocation. Used by the per-conformer
/// planarity acceptance check ([`crate::checks::planarity_ok`]).
pub fn improper_energy(spec: &MoleculeSpec, coords: &[f64], dim: usize) -> f64 {
    let mut e = 0.0;
    for imp in &spec.impropers {
        if let Some(gm) = improper_geom(coords, imp, dim) {
            e += improper_energy_one(imp, &gm);
        }
    }
    e
}

pub fn improper_energy_grad(spec: &MoleculeSpec, coords: &[f64], dim: usize) -> (f64, Vec<f64>) {
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];
    for imp in &spec.impropers {
        let Some(gm) = improper_geom(coords, imp, dim) else {
            continue;
        };
        let (i, j, k, l) = (
            imp.atoms[0] as usize,
            imp.atoms[1] as usize,
            imp.atoms[2] as usize,
            imp.atoms[3] as usize,
        );
        e += improper_energy_one(imp, &gm);
        let ImproperGeom { rji, rjk, rjl, dji, djk, djl, cosy, siny_sq } = gm;
        let (c1, c2, fc) = (imp.c1, imp.c2, imp.fc);

        // gradient — near a singular geometry (L along the plane normal → sinY→0, or collinear
        // i-j-k → sinTheta→0) RDKit does NOT skip: it floors sinY and sinTheta at 1e-8
        // (`std::max(sqrt(...), 1e-8)` in InversionContribs::getGrad) and computes the large but
        // finite gradient. We match that exactly — skipping produced up to rel~1.6 vs RDKit on
        // pathological corpus conformers, and made native's minimizer follow a different gradient
        // than RDKit's there.
        let cos_theta = dot3(rji, rjk).clamp(-1.0, 1.0);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let siny = siny_sq.sqrt().max(1e-8);
        let sin_theta = sin_theta_sq.sqrt().max(1e-8);
        let de_dw = -fc * (c1 * cosy + 4.0 * c2 * cosy * siny);
        let t1 = cross(rjl, rjk);
        let t2 = cross(rji, rjl);
        let t3 = cross(rjk, rji);
        let term1 = siny * sin_theta;
        let term2 = cosy / (siny * sin_theta_sq);
        let mut tg1 = [0.0; 3];
        let mut tg3 = [0.0; 3];
        let mut tg4 = [0.0; 3];
        for c in 0..3 {
            tg1[c] = (t1[c] / term1 - (rji[c] - rjk[c] * cos_theta) * term2) / dji;
            tg3[c] = (t2[c] / term1 - (rjk[c] - rji[c] * cos_theta) * term2) / djk;
            tg4[c] = (t3[c] / term1 - rjl[c] * cosy / siny) / djl;
        }
        for c in 0..3 {
            g[i * dim + c] += de_dw * tg1[c];
            g[j * dim + c] += -de_dw * (tg1[c] + tg3[c] + tg4[c]);
            g[k * dim + c] += de_dw * tg3[c];
            g[l * dim + c] += de_dw * tg4[c];
        }
    }
    (e, g)
}

/// Stage-C distance constraint: flat-bottom harmonic pinning `|p_i − p_j|` into `[min_len, max_len]`.
/// `E = 0.5·fc·(d − bound)²` outside the window, 0 inside (RDKit `DistanceConstraintContrib`).
#[derive(Clone, Debug)]
pub struct DistConstraint {
    pub i: usize,
    pub j: usize,
    pub min_len: f64,
    pub max_len: f64,
    pub fc: f64,
}

const KNOWN_DIST_TOL: f64 = 0.01;
const KNOWN_DIST_FORCE_CONSTANT: f64 = 100.0;

/// Build the Stage-C constraint set from `coords` (dim 3, the post-Stage-A geometry). Each atom pair
/// is constrained at most once, in this order (RDKit `construct3DForceField`):
///
/// - 1-4 pairs carrying an experimental torsion: no distance constraint (the torsion term covers them).
/// - 1-2 pairs (`spec.bonds`): pinned to their current distance ±0.01 Å, force constant 100.
/// - 1-3 pairs (`spec.angles`) marked `triple`: a 179–180° [`AngleConstraint`], force constant 1.
/// - 1-3 pairs centered on an improper center: the bounds-matrix limits, force constant 100.
/// - all other 1-3 pairs: pinned to their current distance ±0.01 Å, force constant 100.
/// - every remaining pair: the bounds-matrix limits, force constant `bounds_force_scaling × 10`.
pub fn build_stage_c_constraints(
    spec: &MoleculeSpec,
    coords: &[f64],
) -> (Vec<DistConstraint>, Vec<AngleConstraint>) {
    let n = spec.n_atoms;
    let mut marked = vec![false; n * n];
    let mut mark = |a: usize, b: usize| {
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        marked[lo * n + hi] = true;
    };
    let dist = |i: usize, j: usize| -> f64 {
        let (pi, pj) = (pt(coords, i, 3), pt(coords, j, 3));
        ((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2)).sqrt()
    };
    let mut cs = Vec::new();

    // torsions only mark their 1-4 pair (the torsion energy handles them)
    for t in &spec.exp_torsions {
        mark(t.atoms[0] as usize, t.atoms[3] as usize);
    }
    // 1-2 (bonds): pin to current distance
    for b in &spec.bonds {
        let (i, j) = (b[0] as usize, b[1] as usize);
        let d = dist(i, j);
        cs.push(DistConstraint {
            i,
            j,
            min_len: d - KNOWN_DIST_TOL,
            max_len: d + KNOWN_DIST_TOL,
            fc: KNOWN_DIST_FORCE_CONSTANT,
        });
        mark(i, j);
    }
    // improper centers → their 1-3 distances come from the bounds matrix (add13Terms
    // `isImproperConstrained` branch), so the inversion term, not a stale distance pin, sets planarity.
    let mut is_imp_center = vec![false; n];
    for imp in &spec.impropers {
        is_imp_center[imp.atoms[1] as usize] = true;
    }
    // 1-3 (angles)
    let mut angle_cs = Vec::new();
    for a in &spec.angles {
        let (i, jc, k) = (
            a.atoms[0] as usize,
            a.atoms[1] as usize,
            a.atoms[2] as usize,
        );
        if a.triple {
            // triple bonds / allenes: a 179-180° angle constraint (not a distance pin)
            angle_cs.push(AngleConstraint {
                i,
                j: jc,
                k,
                min_angle: 179.0,
                max_angle: 180.0,
                fc: 1.0,
            });
        } else if is_imp_center[jc] {
            cs.push(DistConstraint {
                i,
                j: k,
                min_len: spec.lb64(i, k),
                max_len: spec.ub64(i, k),
                fc: KNOWN_DIST_FORCE_CONSTANT,
            });
        } else {
            let d = dist(i, k);
            cs.push(DistConstraint {
                i,
                j: k,
                min_len: d - KNOWN_DIST_TOL,
                max_len: d + KNOWN_DIST_TOL,
                fc: KNOWN_DIST_FORCE_CONSTANT,
            });
        }
        mark(i, k);
    }
    // long-range: every remaining pair constrained to the bounds matrix
    let fc_long = spec.bounds_force_scaling as f64 * 10.0;
    for i in 0..n {
        for j in (i + 1)..n {
            if !marked[i * n + j] {
                cs.push(DistConstraint {
                    i,
                    j,
                    min_len: spec.lb64(i, j),
                    max_len: spec.ub64(i, j),
                    fc: fc_long,
                });
            }
        }
    }
    (cs, angle_cs)
}

/// One flat-bottom distance constraint's contribution (RDKit `DistanceConstraintContribs`:
/// `E = 0.5·fc·diff²`, `grad = fc·(d−bound)/d·(p_i−p_j)`, zero inside `[min,max]`).
#[inline]
fn accum_dist_constraint(c: &DistConstraint, coords: &[f64], e: &mut f64, g: &mut [f64]) {
    let (pi, pj) = (pt(coords, c.i, 3), pt(coords, c.j, 3));
    let d = ((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2)).sqrt();
    let bound = if d < c.min_len {
        c.min_len
    } else if d > c.max_len {
        c.max_len
    } else {
        return;
    };
    let diff = d - bound; // signed
    *e += 0.5 * c.fc * diff * diff;
    let pre = c.fc * diff / d.max(1e-8);
    for t in 0..3 {
        let gg = pre * (pi[t] - pj[t]);
        g[c.i * 3 + t] += gg;
        g[c.j * 3 + t] -= gg;
    }
}

/// Energy + gradient of a Stage-C distance-constraint set (3D coords), RDKit
/// `DistanceConstraintContribs`: `E = Σ 0.5·fc·diff²`, `grad = fc·(d−bound)/d·(p_i−p_j)`.
pub fn dist_constraint_energy_grad(cs: &[DistConstraint], coords: &[f64]) -> (f64, Vec<f64>) {
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];
    for c in cs {
        accum_dist_constraint(c, coords, &mut e, &mut g);
    }
    (e, g)
}

/// Stage-C angle constraint: flat-bottom harmonic on the angle i-j-k (j central), in **degrees**.
/// Built only for near-linear geometries (triple bonds, allenes).
#[derive(Clone, Debug)]
pub struct AngleConstraint {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub min_angle: f64,
    pub max_angle: f64,
    pub fc: f64,
}

/// Energy + gradient of Stage-C angle constraints (3D), RDKit `AngleConstraintContribs`:
/// `E = Σ fc·angleTerm²` with `angleTerm` the degrees outside `[min,max]`.
pub fn angle_constraint_energy_grad(cs: &[AngleConstraint], coords: &[f64]) -> (f64, Vec<f64>) {
    const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];
    for c in cs {
        let (p1, p2, p3) = (pt(coords, c.i, 3), pt(coords, c.j, 3), pt(coords, c.k, 3));
        let r0 = subv(p1, p2);
        let r1 = subv(p3, p2);
        let l0 = dot3(r0, r0).max(1e-5);
        let l1 = dot3(r1, r1).max(1e-5);
        let cos = (dot3(r0, r1) / (l0 * l1).sqrt()).clamp(-1.0, 1.0);
        let angle = RAD2DEG * cos.acos();
        let angle_term = if angle < c.min_angle {
            angle - c.min_angle
        } else if angle > c.max_angle {
            angle - c.max_angle
        } else {
            continue;
        };
        e += c.fc * angle_term * angle_term;
        let de_dtheta = 2.0 * RAD2DEG * c.fc * angle_term;
        let rp = cross(r1, r0);
        let prefactor = de_dtheta / norm3(rp).max(1e-5);
        let (t0, t1) = (-prefactor / l0, prefactor / l1);
        let dedp0 = scale3(cross(r0, rp), t0);
        let dedp2 = scale3(cross(r1, rp), t1);
        for c3 in 0..3 {
            g[c.i * 3 + c3] += dedp0[c3];
            g[c.k * 3 + c3] += dedp2[c3];
            g[c.j * 3 + c3] += -dedp0[c3] - dedp2[c3];
        }
    }
    (e, g)
}

/// Stage-C objective (3D): M6 torsions + UFF impropers + the distance and angle constraints
/// pre-built by [`build_stage_c_constraints`]. No chiral term.
pub fn stage_c_energy_grad(
    spec: &MoleculeSpec,
    dist_c: &[DistConstraint],
    angle_c: &[AngleConstraint],
    coords: &[f64],
) -> (f64, Vec<f64>) {
    ff_tick();
    let (mut e, mut g) = torsion_energy_grad(spec, coords, 3);
    let (ed, gd) = dist_constraint_energy_grad(dist_c, coords);
    let (ei, gi) = improper_energy_grad(spec, coords, 3);
    let (ea, ga) = angle_constraint_energy_grad(angle_c, coords);
    e += ed + ei + ea;
    for c in 0..g.len() {
        g[c] += gd[c] + gi[c] + ga[c];
    }
    (e, g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bb_core::{ChiralSet, ExpTorsion, Improper, MoleculeSpec};

    fn coords(n: usize, dim: usize, mut s: u64) -> Vec<f64> {
        (0..n * dim)
            .map(|_| {
                s = s.wrapping_add(0x9E3779B97F4A7C15);
                let mut z = s;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
                ((z ^ (z >> 31)) as f64 / u64::MAX as f64 - 0.5) * 6.0
            })
            .collect()
    }

    fn synthetic_spec() -> MoleculeSpec {
        let n = 5;
        let mut bounds = vec![0.0f32; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                bounds[i * n + j] = 3.0; // UB
                bounds[j * n + i] = 1.3; // LB
            }
        }
        MoleculeSpec {
            n_atoms: n,
            dim: 4,
            bounds,
            chiral_sets: vec![ChiralSet {
                center: 4,
                atoms: [0, 1, 2, 3],
                vol_lo: 5.0,
                vol_hi: 100.0,
                fused_small_rings: false,
            }],
            ..Default::default()
        }
    }

    // The hot path (minimizer) uses the pre-built-pair evals; the reject/certify path uses the basin
    // scan. They MUST agree bit-for-bit or the ensembles the minimizer produces would diverge from the
    // RDKit-certified gradient. Assert byte-identity of energy and every gradient component across both
    // basin regimes (BASIN_ALL keeps all pairs; BASIN_DEFAULT exercises the filter), for Stage A and B.
    #[test]
    fn pair_eval_bit_identical_to_scan() {
        let spec = synthetic_spec();
        let dim = 4;
        let x = coords(spec.n_atoms, dim, 0xC0FFEE);
        for basin in [BASIN_ALL, BASIN_DEFAULT] {
            let pairs = build_dist_pairs(&spec, basin);
            let (ea, ga) = stage_a_energy_grad(&spec, &x, dim, basin);
            let (ep, gp) = stage_a_energy_grad_pairs(&spec, &x, dim, &pairs);
            assert_eq!(ea.to_bits(), ep.to_bits(), "Stage-A energy differs (basin {basin})");
            assert_eq!(ga, gp, "Stage-A gradient differs (basin {basin})");
            let (eb, gb) = stage_b_energy_grad(&spec, &x, dim, basin);
            let (ebp, gbp) = stage_b_energy_grad_pairs(&spec, &x, dim, &pairs);
            assert_eq!(eb.to_bits(), ebp.to_bits(), "Stage-B energy differs (basin {basin})");
            assert_eq!(gb, gbp, "Stage-B gradient differs (basin {basin})");
        }
    }

    #[test]
    fn gradient_matches_finite_difference() {
        let spec = synthetic_spec();
        let dim = 4;
        let x = coords(spec.n_atoms, dim, 0xC0FFEE);
        let (_e, g) = stage_a_energy_grad(&spec, &x, dim, BASIN_ALL);

        let eps = 1e-6;
        let mut max_rel = 0.0f64;
        for idx in 0..x.len() {
            let mut xp = x.clone();
            let mut xm = x.clone();
            xp[idx] += eps;
            xm[idx] -= eps;
            let fd = (stage_a_energy_grad(&spec, &xp, dim, BASIN_ALL).0
                - stage_a_energy_grad(&spec, &xm, dim, BASIN_ALL).0)
                / (2.0 * eps);
            let denom = g[idx].abs().max(fd.abs()).max(1e-6);
            max_rel = max_rel.max((g[idx] - fd).abs() / denom);
        }
        assert!(
            max_rel < 1e-4,
            "analytic vs finite-diff: max rel err {max_rel:.2e}"
        );
    }

    #[test]
    fn torsion_gradient_matches_finite_difference() {
        let spec = MoleculeSpec {
            n_atoms: 4,
            dim: 3,
            bounds: vec![0.0f32; 16],
            chiral_sets: vec![],
            tetrahedral_centers: vec![],
            // Only the m=1..4 force constants are exercised here: RDKit's analytic gradient is
            // deliberately inconsistent with its own energy for the higher terms — the m=6 term
            // reuses forceConstants[4]/signs[4] (not [5]) and drops a cosφ factor (6·sinφ, not
            // 6·cosφ·sinφ), and there is no gradient term reading forceConstants[5] at all. We
            // replicate that quirk for bit-exact parity, so `gradient == d(energy)/dx` can only hold
            // where V[4]=V[5]=0. The full six-term gradient (quirk included) is validated directly
            // against RDKit's TorsionAngleContribs in bb-rdkit/tests/ff_parity.rs (rel < 1e-5).
            exp_torsions: vec![ExpTorsion {
                atoms: [0, 1, 2, 3],
                v: [1.0, 0.5, 0.3, 0.2, 0.0, 0.0],
                signs: [1, -1, 1, -1, 1, -1],
            }],
            ..Default::default()
        };
        let dim = 3;
        // three random geometries so we don't accidentally sit at a sin(mφ)=0 special point
        for seed in [0xBEEFu64, 0x1234, 0xF00D] {
            let x = coords(spec.n_atoms, dim, seed);
            let (_e, g) = torsion_energy_grad(&spec, &x, dim);
            let eps = 1e-6;
            let mut max_rel = 0.0f64;
            for idx in 0..x.len() {
                let mut xp = x.clone();
                let mut xm = x.clone();
                xp[idx] += eps;
                xm[idx] -= eps;
                let fd = (torsion_energy_grad(&spec, &xp, dim).0
                    - torsion_energy_grad(&spec, &xm, dim).0)
                    / (2.0 * eps);
                let denom = g[idx].abs().max(fd.abs()).max(1e-6);
                max_rel = max_rel.max((g[idx] - fd).abs() / denom);
            }
            assert!(
                max_rel < 1e-4,
                "seed {seed:#x}: torsion grad vs fd max rel err {max_rel:.2e}"
            );
        }
    }

    #[test]
    fn dist_constraint_gradient_matches_finite_difference() {
        let dim = 3;
        let x = coords(5, dim, 0xABCD);
        // mix of under-min (active), over-max (active), and inside (inactive) constraints
        let cs = vec![
            DistConstraint {
                i: 0,
                j: 1,
                min_len: 5.0,
                max_len: 6.0,
                fc: 100.0,
            },
            DistConstraint {
                i: 1,
                j: 2,
                min_len: 0.05,
                max_len: 0.5,
                fc: 50.0,
            },
            DistConstraint {
                i: 2,
                j: 3,
                min_len: 0.0,
                max_len: 100.0,
                fc: 10.0,
            },
            DistConstraint {
                i: 3,
                j: 4,
                min_len: 1.0,
                max_len: 2.0,
                fc: 20.0,
            },
        ];
        let (_e, g) = dist_constraint_energy_grad(&cs, &x);
        let eps = 1e-6;
        let mut max_rel = 0.0f64;
        for idx in 0..x.len() {
            let mut xp = x.clone();
            let mut xm = x.clone();
            xp[idx] += eps;
            xm[idx] -= eps;
            let fd = (dist_constraint_energy_grad(&cs, &xp).0
                - dist_constraint_energy_grad(&cs, &xm).0)
                / (2.0 * eps);
            let denom = g[idx].abs().max(fd.abs()).max(1e-6);
            max_rel = max_rel.max((g[idx] - fd).abs() / denom);
        }
        assert!(
            max_rel < 1e-4,
            "dist-constraint grad vs fd max rel err {max_rel:.2e}"
        );
    }

    #[test]
    fn angle_constraint_gradient_matches_finite_difference() {
        let dim = 3;
        // random 3-atom geometry; the angle is generically < 179°, so the [179,180] constraint is active
        for seed in [0x5u64, 0x6, 0x7] {
            let x = coords(3, dim, seed);
            let cs = vec![AngleConstraint {
                i: 0,
                j: 1,
                k: 2,
                min_angle: 179.0,
                max_angle: 180.0,
                fc: 1.0,
            }];
            let (_e, g) = angle_constraint_energy_grad(&cs, &x);
            let eps = 1e-7;
            let mut max_rel = 0.0f64;
            for idx in 0..x.len() {
                let mut xp = x.clone();
                let mut xm = x.clone();
                xp[idx] += eps;
                xm[idx] -= eps;
                let fd = (angle_constraint_energy_grad(&cs, &xp).0
                    - angle_constraint_energy_grad(&cs, &xm).0)
                    / (2.0 * eps);
                let denom = g[idx].abs().max(fd.abs()).max(1e-6);
                max_rel = max_rel.max((g[idx] - fd).abs() / denom);
            }
            assert!(
                max_rel < 1e-3,
                "seed {seed:#x}: angle-constraint grad vs fd max rel err {max_rel:.2e}"
            );
        }
    }

    #[test]
    fn improper_gradient_matches_finite_difference() {
        let spec = MoleculeSpec {
            n_atoms: 4,
            dim: 3,
            bounds: vec![0.0f32; 16],
            impropers: vec![Improper {
                atoms: [0, 1, 2, 3], // center = atoms[1]
                c0: 1.0,
                c1: -1.0,
                c2: 0.5,
                fc: 30.0,
            }],
            ..Default::default()
        };
        let dim = 3;
        for seed in [0x11u64, 0x22, 0x33, 0x44] {
            let x = coords(4, dim, seed);
            let (_e, g) = improper_energy_grad(&spec, &x, dim);
            let eps = 1e-7;
            let mut max_rel = 0.0f64;
            for idx in 0..x.len() {
                let mut xp = x.clone();
                let mut xm = x.clone();
                xp[idx] += eps;
                xm[idx] -= eps;
                let fd = (improper_energy_grad(&spec, &xp, dim).0
                    - improper_energy_grad(&spec, &xm, dim).0)
                    / (2.0 * eps);
                let denom = g[idx].abs().max(fd.abs()).max(1e-6);
                max_rel = max_rel.max((g[idx] - fd).abs() / denom);
            }
            assert!(
                max_rel < 1e-3,
                "seed {seed:#x}: improper grad vs fd max rel err {max_rel:.2e}"
            );
        }
    }
}
