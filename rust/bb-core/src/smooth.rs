//! Triangle-inequality smoothing of a distance-bounds matrix — the one distance-geometry primitive
//! shared by every consumer of [`crate::MoleculeSpec::bounds`]: the RDKit-free spec builder
//! (`bb-spec`) smooths the freshly built matrix, and the embed engine (`bb-embed`) re-smooths its
//! working matrix after pinning. It lives here, next to the matrix it operates on, so neither of those
//! crates has to depend on the other for it.
//!
//! Storage matches [`crate::MoleculeSpec::bounds`]: `n×n` row-major, `b[r*n+c]` with `r<c` = upper
//! bound, `r>c` = lower bound.

/// Triangle-inequality bounds smoothing in full `f64`, cast once by the caller (RDKit
/// `triangleSmoothBounds`, which operates on the `double` bounds matrix). Structurally identical to
/// [`triangle_smooth`] but with no intermediate rounding: the bridge (`bridge.cc`) smooths the
/// `double` matrix and only then casts to `f32` for the spec, so the RDKit-free spec path must do the
/// same to reproduce those `f32` bounds bit-for-bit. Smoothing the already-cast `f32` matrix instead
/// leaves ~1 ULP differences in a fraction of cells, which would perturb the stochastic embed.
pub fn triangle_smooth_f64(b: &mut [f64], n: usize, tol: f64) -> bool {
    if n < 2 {
        return true;
    }
    for k in 0..n {
        for i in 0..(n - 1) {
            if i == k {
                continue;
            }
            let (ii, ik) = if i < k { (i, k) } else { (k, i) };
            let uik = b[ii * n + ik]; // upper(i,k)
            let lik = b[ik * n + ii]; // lower(i,k)
            for j in (i + 1)..n {
                if j == k {
                    continue;
                }
                let (jj, jk) = if j < k { (j, k) } else { (k, j) };
                let ukj = b[jj * n + jk]; // upper(k,j)

                // tighten upper(i,j) = b[i*n+j]
                let sum = uik + ukj;
                if b[i * n + j] > sum {
                    b[i * n + j] = sum;
                }

                // raise lower(i,j) = b[j*n+i] (exact RDKit if/else on the ORIGINAL value)
                let diff_lik_ukj = lik - ukj;
                let diff_ljk_uik = b[jk * n + jj] - uik; // lower(k,j) − upper(i,k)
                let cur_lower = b[j * n + i];
                if cur_lower < diff_lik_ukj {
                    b[j * n + i] = diff_lik_ukj;
                } else if cur_lower < diff_ljk_uik {
                    b[j * n + i] = diff_ljk_uik;
                }

                let l_bound = b[j * n + i];
                let u_bound = b[i * n + j];
                if tol > 0.0
                    && (l_bound - u_bound) / l_bound > 0.0
                    && (l_bound - u_bound) / l_bound < tol
                {
                    b[i * n + j] = l_bound; // within tol: bump UB up to LB
                } else if l_bound - u_bound > 0.0 {
                    return false; // inconsistent
                }
            }
        }
    }
    true
}

/// Triangle-inequality bounds smoothing (RDKit `triangleSmoothBounds`): a Floyd–Warshall pass over
/// `k` tightening every upper bound and raising every lower bound. When a lower bound ends up above
/// its upper bound by less than the fraction `tol`, the upper bound is raised to meet it. Returns
/// false if a pair's bounds are inconsistent by more than `tol`.
///
/// This `f32` form is what the embed engine runs on its working matrix; for spec construction, where
/// the goal is to reproduce the bridge's `f32` bounds bit-for-bit, use [`triangle_smooth_f64`] and
/// cast once.
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
}
