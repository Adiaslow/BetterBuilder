//! DistGeom bounds-matrix operations: coordMap tightening + triangle-inequality smoothing — both
//! standard distance-geometry algorithm logic. Used by the two-stage recipe: for each core seed the
//! frozen-core pair distances are written into the bounds and propagated to every atom by
//! re-smoothing, matching RDKit's `adjustBoundsMatFromCoordMap` + `triangleSmoothBounds`.
//!
//! Storage matches [`bb_core::MoleculeSpec::bounds`]: `n×n` row-major, `b[r*n+c]` with `r<c` = upper
//! bound, `r>c` = lower bound.

/// Set lower = upper = |coord_i − coord_j| for every pair of pinned atoms
/// (RDKit `adjustBoundsMatFromCoordMap`). `pinned[a] = Some(xyz)` marks a pinned atom.
pub fn adjust_from_coord_map(bounds: &mut [f32], n: usize, pinned: &[Option<[f64; 3]>]) {
    // Carry each pinned atom's coordinate alongside its index, so the pair loop needs no unwrap.
    let idx: Vec<(usize, [f64; 3])> = (0..pinned.len())
        .filter_map(|a| pinned[a].map(|p| (a, p)))
        .collect();
    for a in 0..idx.len() {
        for b in (a + 1)..idx.len() {
            let (i, pi) = idx[a];
            let (j, pj) = idx[b];
            let d = ((pi[0] - pj[0]).powi(2) + (pi[1] - pj[1]).powi(2) + (pi[2] - pj[2]).powi(2))
                .sqrt();
            let (lo, hi) = if i < j { (i, j) } else { (j, i) };
            bounds[lo * n + hi] = d as f32; // upper
            bounds[hi * n + lo] = d as f32; // lower
        }
    }
}

/// Triangle-inequality bounds smoothing — an exact port of RDKit `triangleSmoothBounds` (Floyd–
/// Warshall over `k`; note the lower-bound `if/else` and the `tol` rule that bumps UB up to LB when
/// they cross by less than `tol` fraction). Returns false if the bounds become inconsistent.
pub fn triangle_smooth(b: &mut [f32], n: usize, tol: f64) -> bool {
    if n < 2 {
        return true;
    }
    for k in 0..n {
        for i in 0..(n - 1) {
            if i == k {
                continue;
            }
            let (ii, ik) = if i < k { (i, k) } else { (k, i) };
            let uik = b[ii * n + ik] as f64; // upper(i,k)
            let lik = b[ik * n + ii] as f64; // lower(i,k)
            for j in (i + 1)..n {
                if j == k {
                    continue;
                }
                let (jj, jk) = if j < k { (j, k) } else { (k, j) };
                let ukj = b[jj * n + jk] as f64; // upper(k,j)

                // tighten upper(i,j) = b[i*n+j]
                let sum = uik + ukj;
                if (b[i * n + j] as f64) > sum {
                    b[i * n + j] = sum as f32;
                }

                // raise lower(i,j) = b[j*n+i] (exact RDKit if/else on the ORIGINAL value)
                let diff_lik_ukj = lik - ukj;
                let diff_ljk_uik = (b[jk * n + jj] as f64) - uik; // lower(k,j) − upper(i,k)
                let cur_lower = b[j * n + i] as f64;
                if cur_lower < diff_lik_ukj {
                    b[j * n + i] = diff_lik_ukj as f32;
                } else if cur_lower < diff_ljk_uik {
                    b[j * n + i] = diff_ljk_uik as f32;
                }

                let l_bound = b[j * n + i] as f64;
                let u_bound = b[i * n + j] as f64;
                if tol > 0.0
                    && (l_bound - u_bound) / l_bound > 0.0
                    && (l_bound - u_bound) / l_bound < tol
                {
                    b[i * n + j] = l_bound as f32; // within tol: bump UB up to LB
                } else if l_bound - u_bound > 0.0 {
                    return false; // inconsistent
                }
            }
        }
    }
    true
}

/// Sidechain bounds for one core seed: clone `base`, pin the frozen-core distances, re-smooth
/// (tol 0.05). On smoothing failure, fall back to the base bounds (rely on the freeze; RDKit instead
/// rebuilds without 1-5 bounds — a rare path, since a seed that passed its checks is bounds-consistent).
pub fn coord_map_bounds(base: &[f32], n: usize, pinned: &[Option<[f64; 3]>]) -> Vec<f32> {
    let mut b = base.to_vec();
    adjust_from_coord_map(&mut b, n, pinned);
    if triangle_smooth(&mut b, n, 0.05) {
        b
    } else {
        base.to_vec()
    }
}

#[cfg(test)]
mod tests {
    // Tests index the n×n bounds matrix in explicit `row * n + col` form to document its layout
    // (e.g. `b[0 * n + 2]` = upper(0,2)); the `0 *`/`+ 0` terms are intentional, not dead arithmetic.
    #![allow(clippy::identity_op, clippy::erasing_op)]
    use super::*;

    #[test]
    fn smoothing_tightens_upper_via_triangle() {
        // 3 atoms: UB(0,1)=1, UB(1,2)=1, UB(0,2)=5 -> should tighten to 2.
        let n = 3;
        let mut b = vec![0.0f32; n * n];
        // upper (r<c)
        b[0 * n + 1] = 1.0;
        b[1 * n + 2] = 1.0;
        b[0 * n + 2] = 5.0;
        // lower (r>c)
        b[1 * n + 0] = 0.5;
        b[2 * n + 1] = 0.5;
        b[2 * n + 0] = 0.5;
        assert!(triangle_smooth(&mut b, n, 0.0));
        assert!(
            (b[0 * n + 2] - 2.0).abs() < 1e-5,
            "UB(0,2) not tightened: {}",
            b[0 * n + 2]
        );
    }

    #[test]
    fn coord_map_pins_distances() {
        let n = 3;
        let mut b = vec![0.0f32; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                b[i * n + j] = 10.0; // loose upper
                b[j * n + i] = 0.1; // loose lower
            }
        }
        let pinned = vec![Some([0.0, 0.0, 0.0]), Some([3.0, 0.0, 0.0]), None];
        adjust_from_coord_map(&mut b, n, &pinned);
        assert!((b[0 * n + 1] - 3.0).abs() < 1e-5); // upper(0,1) pinned to 3
        assert!((b[1 * n + 0] - 3.0).abs() < 1e-5); // lower(0,1) pinned to 3
    }
}
