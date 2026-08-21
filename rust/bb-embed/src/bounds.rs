//! DistGeom bounds-matrix operations specific to the core-pin embed: coordMap tightening
//! ([`adjust_from_coord_map_f64`]) combined with re-smoothing by [`coord_map_bounds_f64`]. The
//! triangle-inequality primitive itself lives in [`bb_core::smooth`], shared with the spec builder.
//!
//! Storage matches [`bb_core::MoleculeSpec::bounds`]: `n×n` row-major, `b[r*n+c]` with `r<c` = upper
//! bound, `r>c` = lower bound.

use bb_core::smooth::triangle_smooth_f64;

/// Triangle-smoothing tolerance for the coordMap re-smooth — RDKit's `DGeomHelpers` default (`tol`
/// argument to `triangleSmoothBounds`); also the `setupInitialBoundsMatrix` smoothing tol.
const SMOOTH_TOL: f64 = 0.05;

/// Set lower = upper = |coord_i − coord_j| for every pair of pinned atoms
/// (RDKit `adjustBoundsMatFromCoordMap`). `pinned[a] = Some(xyz)` marks a pinned atom. f64 throughout,
/// matching RDKit's coordMap path (the embed reads these via `ub64`; there is no f32 twin — the spec's
/// f32 `bounds` shadow is derived from the tightened f64 in `MoleculeSpec::with_seed_bounds`).
pub fn adjust_from_coord_map_f64(bounds: &mut [f64], n: usize, pinned: &[Option<[f64; 3]>]) {
    let idx: Vec<(usize, [f64; 3])> = (0..pinned.len())
        .filter_map(|a| pinned[a].map(|p| (a, p)))
        .collect();
    for a in 0..idx.len() {
        for b in (a + 1)..idx.len() {
            let (i, pi) = idx[a];
            let (j, pj) = idx[b];
            let d = ((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2)).sqrt();
            let (lo, hi) = if i < j { (i, j) } else { (j, i) };
            bounds[lo * n + hi] = d;
            bounds[hi * n + lo] = d;
        }
    }
}

/// The sidechain bounds RDKit's coordMap path produces, at machine precision (f64 raw bounds, f64
/// pins, f64 smoothing): `base` with the pinned pair distances written in ([`adjust_from_coord_map_f64`])
/// and re-smoothed at [`SMOOTH_TOL`]. Falls back to the pin-free smoothed bounds if the pinned
/// distances are triangle-inconsistent (unreachable for pins taken from a real embedded conformer,
/// whose distances are realizable — the only way the recipe produces them).
pub fn coord_map_bounds_f64(base: &[f64], n: usize, pinned: &[Option<[f64; 3]>]) -> Vec<f64> {
    let mut b = base.to_vec();
    adjust_from_coord_map_f64(&mut b, n, pinned);
    if triangle_smooth_f64(&mut b, n, SMOOTH_TOL) {
        return b;
    }
    let mut b = base.to_vec();
    triangle_smooth_f64(&mut b, n, SMOOTH_TOL);
    b
}

#[cfg(test)]
mod tests {
    // Tests index the n×n bounds matrix in explicit `row * n + col` form to document its layout
    // (e.g. `b[0 * n + 2]` = upper(0,2)); the `0 *`/`+ 0` terms are intentional, not dead arithmetic.
    #![allow(clippy::identity_op, clippy::erasing_op)]
    use super::*;

    #[test]
    fn coord_map_pins_distances() {
        let n = 3;
        let mut b = vec![0.0f64; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                b[i * n + j] = 10.0; // loose upper
                b[j * n + i] = 0.1; // loose lower
            }
        }
        let pinned = vec![Some([0.0, 0.0, 0.0]), Some([3.0, 0.0, 0.0]), None];
        adjust_from_coord_map_f64(&mut b, n, &pinned);
        assert!((b[0 * n + 1] - 3.0).abs() < 1e-9); // upper(0,1) pinned to 3
        assert!((b[1 * n + 0] - 3.0).abs() < 1e-9); // lower(0,1) pinned to 3
    }
}
