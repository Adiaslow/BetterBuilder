//! Stage-A DistGeom force field: distance-bound violation (all dims) + chiral-volume (x,y,z) +
//! 4th-dimension penalty. f64, analytic gradients, finite-difference checked. This is the
//! candidate's reimplementation of ETKDG's Stage-A objective; ultimate validation is ensemble
//! RMSD vs the oracle. To stay faithful, we **match RDKit's actual gradients** — RDKit's
//! `ChiralViolationContribs`/`FourthDimContribs` omit the factor of 2 (their gradient is the
//! derivative of `0.5·w·x²`, not `w·x²`), so we use energy `0.5·w·x²` + gradient `w·x`, giving the
//! exact same force RDKit descends on (and keeping energy/gradient consistent for finite-diff).
//! The distance term already copies RDKit's exact preFactors.

use std::sync::atomic::{AtomicU64, Ordering};

use bb_core::MoleculeSpec;

/// Diagnostic counter: total DistGeom (Stage A/B) + Stage-C force-field evaluations. Gated behind the
/// `profile` feature so it adds zero overhead (and no cross-thread atomic contention) in production /
/// parallel runs; enable with `--features profile` for profiling.
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

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
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

/// Distance-bound violation term (all dims), accumulated into `e`/`g` — the O(N²) DistGeom pairwise
/// energy. RDKit's exact 4/8 preFactors; the basin filter drops loosely-bounded (flexible) pairs.
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
            let ub = spec.ub(i, j) as f64;
            let lb = spec.lb(i, j) as f64;
            if ub - lb > basin {
                continue; // basin filter: skip loosely-bounded (flexible) pairs
            }
            let (ub2, lb2) = (ub * ub, lb * lb);
            let mut d2 = 0.0;
            for c in 0..dim {
                let d = coords[i * dim + c] - coords[j * dim + c];
                d2 += d * d;
            }
            let (val, pre) = if d2 > ub2 {
                (d2 / ub2 - 1.0, 4.0 * (d2 / ub2 - 1.0) * (d2.sqrt() / ub2))
            } else if d2 < lb2 {
                let s = d2 + lb2;
                (
                    2.0 * lb2 / s - 1.0,
                    8.0 * lb2 * d2.sqrt() * (1.0 - 2.0 * lb2 / s) / (s * s),
                )
            } else {
                (0.0, 0.0)
            };
            if val > 0.0 {
                *e += W_DIST * val * val;
                let d = d2.sqrt();
                if d > 1e-12 {
                    for c in 0..dim {
                        let gg = W_DIST * pre * (coords[i * dim + c] - coords[j * dim + c]) / d;
                        g[i * dim + c] += gg;
                        g[j * dim + c] -= gg;
                    }
                }
            }
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

    // --- chiral-volume violation (x,y,z only), weight 1.0 ---
    for cs in &spec.chiral_sets {
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
        let bound = if vol < cs.vol_lo as f64 {
            Some(cs.vol_lo as f64)
        } else if vol > cs.vol_hi as f64 {
            Some(cs.vol_hi as f64)
        } else {
            None
        };
        if let Some(b) = bound {
            e += 0.5 * w_chiral * (vol - b) * (vol - b);
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

    // --- fourth-dimension penalty (4th coord only), weight 0.1 ---
    if dim == 4 {
        for i in 0..n {
            let w = coords[i * dim + 3];
            e += 0.5 * w_fourth * w * w;
            g[i * dim + 3] += w_fourth * w; // matches RDKit FourthDimContribs (no factor of 2)
        }
    }

    (e, g)
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
fn subv(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
fn pt(p: &[f64], i: usize, dim: usize) -> [f64; 3] {
    [p[i * dim], p[i * dim + 1], p[i * dim + 2]]
}

/// M6 experimental-torsion term (x,y,z only): `E = Σ_{m=1..6} V[m-1]·(1 + signs[m-1]·cos(m·φ))`,
/// φ = signed dihedral(i,j,k,l). Analytic gradient via the standard four-atom dihedral projection
/// (translation-invariant: the four position gradients sum to zero). Carries Divya's patched V.
pub fn torsion_energy_grad(spec: &MoleculeSpec, coords: &[f64], dim: usize) -> (f64, Vec<f64>) {
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];
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
        let b1 = [pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]];
        let b2 = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
        let b3 = [pl[0] - pk[0], pl[1] - pk[1], pl[2] - pk[2]];
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let (n1sq, n2sq, b2sq) = (dot3(n1, n1), dot3(n2, n2), dot3(b2, b2));
        if n1sq < 1e-10 || n2sq < 1e-10 || b2sq < 1e-10 {
            continue; // degenerate (collinear) — skip
        }
        let b2n = b2sq.sqrt();
        let phi = (dot3(cross(n1, n2), b2) / b2n).atan2(dot3(n1, n2));

        let mut de_dphi = 0.0;
        for idx in 0..6 {
            let m = (idx + 1) as f64;
            let v = t.v[idx] as f64;
            let s = t.signs[idx] as f64;
            e += v * (1.0 + s * (m * phi).cos());
            de_dphi += v * s * (-m) * (m * phi).sin();
        }

        let dphi_di = scale3(n1, -b2n / n1sq);
        let dphi_dl = scale3(n2, b2n / n2sq);
        let pp = dot3(b1, b2) / b2sq;
        let qq = dot3(b3, b2) / b2sq;
        let dphi_dj = [
            -(pp + 1.0) * dphi_di[0] + qq * dphi_dl[0],
            -(pp + 1.0) * dphi_di[1] + qq * dphi_dl[1],
            -(pp + 1.0) * dphi_di[2] + qq * dphi_dl[2],
        ];
        let dphi_dk = [
            pp * dphi_di[0] - (qq + 1.0) * dphi_dl[0],
            pp * dphi_di[1] - (qq + 1.0) * dphi_dl[1],
            pp * dphi_di[2] - (qq + 1.0) * dphi_dl[2],
        ];
        for c in 0..3 {
            g[i * dim + c] += de_dphi * dphi_di[c];
            g[j * dim + c] += de_dphi * dphi_dj[c];
            g[k * dim + c] += de_dphi * dphi_dk[c];
            g[l * dim + c] += de_dphi * dphi_dl[c];
        }
    }
    (e, g)
}

/// UFF out-of-plane / inversion term — the ETKDG Stage-C planarity contribution (`InversionContribs`,
/// parameterized by RDKit's UFF inversion coefficients). Per contrib on atoms (i, j=center, k, l):
/// `E = fc·(C0 + C1·sinY + C2·cos2W)`, where `sinY` is the sine of the Wilson out-of-plane angle
/// (`cosY = n̂·r̂JL`, `n̂ ⟂` the i-j-k plane) and `cos2W = 2·sinY² − 1`. Exact RDKit energy + analytic
/// gradient. Coeffs (C0,C1,C2,fc — fc already ×10) come from the spec.
pub fn improper_energy_grad(spec: &MoleculeSpec, coords: &[f64], dim: usize) -> (f64, Vec<f64>) {
    let mut e = 0.0;
    let mut g = vec![0.0f64; coords.len()];
    for imp in &spec.impropers {
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
            continue;
        }
        rji = scale3(rji, 1.0 / dji);
        rjk = scale3(rjk, 1.0 / djk);
        rjl = scale3(rjl, 1.0 / djl);
        // plane normal n̂ = (-r̂JI) × r̂JK
        let mut n = cross(scale3(rji, -1.0), rjk);
        let nn = norm3(n);
        if nn < 1e-8 {
            continue;
        }
        n = scale3(n, 1.0 / nn);
        let cosy = dot3(n, rjl).clamp(-1.0, 1.0);
        let siny_sq = 1.0 - cosy * cosy;
        let (c0, c1, c2, fc) = (imp.c0 as f64, imp.c1 as f64, imp.c2 as f64, imp.fc as f64);

        // energy (RDKit: sinY = sqrt(sinYSq) if >0 else 0)
        let siny_e = if siny_sq > 0.0 { siny_sq.sqrt() } else { 0.0 };
        let cos2w = 2.0 * siny_e * siny_e - 1.0;
        e += fc * (c0 + c1 * siny_e + c2 * cos2w);

        // gradient — skip near-degenerate geometries (collinear i-j-k → sinTheta²→0, or L along the
        // plane normal → sinY→0), where the out-of-plane angle's gradient is singular and would blow
        // up (to inf/NaN) and thrash the line search. Energy stays continuous; other terms move the
        // atom off the singularity. (RDKit's own optimizer caps step size instead; equivalent under
        // the loose ensemble bar.)
        let cos_theta = dot3(rji, rjk).clamp(-1.0, 1.0);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        if siny_sq < 1e-3 || sin_theta_sq < 1e-3 {
            continue;
        }
        let siny = siny_sq.sqrt();
        let sin_theta = sin_theta_sq.sqrt();
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

/// Build the Stage-C distance-constraint set, faithfully replaying RDKit `construct3DForceField`'s
/// `add12Terms` + `add13Terms` + `addLongRangeDistanceConstraints`. Built **once** from the
/// post-Stage-A `coords` (dim=3): 1-2 and 1-3 distances are pinned to the *current* geometry ±tol,
/// everything else to the bounds matrix. `atomPairs` bookkeeping mirrors RDKit so no pair is
/// double-constrained. (Impropers/`isImproperConstrained` and triple-bond angle constraints are not
/// yet modelled — those 1-3s fall back to the current-distance pin, a small deviation.)
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
                min_len: spec.lb(i, k) as f64,
                max_len: spec.ub(i, k) as f64,
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
                    min_len: spec.lb(i, j) as f64,
                    max_len: spec.ub(i, j) as f64,
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

/// Energy + gradient of a Stage-C distance-constraint set (3D coords). Matches RDKit
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
/// RDKit uses these only for near-linear geometries (triple bonds / allenes: 179–180°, fc 1).
#[derive(Clone, Debug)]
pub struct AngleConstraint {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub min_angle: f64,
    pub max_angle: f64,
    pub fc: f64,
}

/// Energy + gradient of Stage-C angle constraints (3D). Exact RDKit `AngleConstraintContribs`:
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

/// Stage-C objective (3D): M6 torsions + UFF impropers + distance constraints (pre-built via
/// [`build_stage_c_constraints`]). No chiral term (RDKit's Stage C has none). This is the full set
/// of terms in RDKit's `construct3DForceField`.
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
            exp_torsions: vec![ExpTorsion {
                atoms: [0, 1, 2, 3],
                v: [1.0, 0.5, 0.3, 0.2, 0.1, 0.05],
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
