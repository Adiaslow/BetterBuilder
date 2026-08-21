//! Initial coordinates for the DistGeom embed:
//! - [`random_box`] — uniform random coordinates (RDKit `useRandomCoords=True`).
//! - [`metric_matrix`] — metric-matrix eigenvalue embedding (RDKit `useRandomCoords=False`).

use bb_core::MoleculeSpec;
use nalgebra::DMatrix;
use rand::Rng;

/// Random-box initialization (ETKDG `useRandomCoords`): each coord uniform in `[-box/2, +box/2]`.
/// (RDKit uses `boxSize=10`, `boxSizeMult=2`.)
pub fn random_box(n_atoms: usize, dim: usize, box_size: f64, rng: &mut impl Rng) -> Vec<f64> {
    (0..n_atoms * dim)
        .map(|_| box_size * (rng.gen::<f64>() - 0.5))
        .collect()
}

const EIGVAL_TOL: f64 = 0.001;

/// Metric (Gram) matrix `T_ij = ½(d_i0² + d_j0² − d_ij²)` from a squared-distance matrix `sq`
/// (`n*n`, symmetric, zero diagonal), with `d_i0²` = squared distance of atom `i` from the centroid.
/// Returns `(T as n*n row-major, d0)`; `d0[i]` is the squared centroid distance `sqD0i` that
/// `computeInitialCoords`' first reject test needs, returned so the caller reuses it rather than
/// recomputing the O(n²) row sums. Steps 2–4 of the metric-matrix method, factored out so the
/// eigenstep can be exercised on a caller-supplied `sq`. (Native embeds via nalgebra's exact
/// `symmetric_eigen`, not RDKit's power iteration, so RDKit's `sumSqD2` power-seed value is not needed.)
pub fn metric_gram(sq: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    // sumSqD2 = (Σ_{i>j} d_ij²) / N²  (lower triangle only, matching RDKit's SymmMatrix data)
    let mut sum_lower = 0.0;
    for i in 0..n {
        for j in 0..i {
            sum_lower += sq[i * n + j];
        }
    }
    let sum_sq_d2 = sum_lower / (n * n) as f64;

    // squared distances from the centroid: d_i0² = (Σ_j d_ij²)/N − sumSqD2
    let mut d0 = vec![0.0f64; n];
    for i in 0..n {
        let row: f64 = (0..n).map(|j| sq[i * n + j]).sum();
        d0[i] = row / n as f64 - sum_sq_d2;
    }

    // metric (Gram) matrix T_ij = ½(d_i0² + d_j0² − d_ij²)
    let mut t = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            t[i * n + j] = 0.5 * (d0[i] + d0[j] - sq[i * n + j]);
        }
    }
    (t, d0)
}

/// The exact metric-matrix embedding from a **given** squared-distance matrix `sq` (`n*n`, symmetric,
/// zero diagonal): steps 2–6 of [`metric_matrix`], factored out so a caller can embed a specific `sq`.
/// Eigendecomposes `T`, takes the top-`dim` eigenpairs by algebraic value, and places `√λ·v`
/// (positive), `0` (~zero), or a random `[-1,1]` (negative) coordinate.
///
/// Returns `None` on RDKit's `computeInitialCoords` reject (for `n > 3`): any squared-centroid
/// distance `sqD0i < EIGVAL_TOL`, or any of the top-`dim` eigenvalues within `EIGVAL_TOL` of zero
/// (`numZeroFail = 1`). The caller re-draws a fresh distance matrix, exactly as RDKit's embed loop
/// does when `computeInitialCoords` returns false. `Some(n*dim coords)` otherwise.
pub fn coords_from_sq(sq: &[f64], n: usize, dim: usize, rng: &mut impl Rng) -> Option<Vec<f64>> {
    // 2–4. metric (Gram) matrix from the squared distances (d0 = sqD0i, reused by the reject below)
    let (t, d0) = metric_gram(sq, n);

    // computeInitialCoords reject 1 (n>3): sqD0i[i] < EIGVAL_TOL.
    if n > 3 && d0.iter().any(|&d| d < EIGVAL_TOL) {
        return None;
    }

    let t = DMatrix::<f64>::from_row_slice(n, n, &t);

    // 5. symmetric eigendecomposition; take the top `dim` eigenpairs (largest eigenvalue first)
    let eig = t.symmetric_eigen();
    let mut order: Vec<usize> = (0..n).collect();
    // Descending eigenvalue order. `total_cmp` gives a total order with no panic path (the eigenvalues
    // of a symmetric matrix are real, so this never sees a NaN, but total_cmp needs no such assumption).
    order.sort_by(|&a, &b| eig.eigenvalues[b].total_cmp(&eig.eigenvalues[a]));

    // computeInitialCoords reject 2 (n>3): reject if any top-`dim` eigenvalue is ~zero (numZeroFail=1).
    if n > 3
        && (0..dim.min(n)).any(|a| eig.eigenvalues[order[a]].abs() < EIGVAL_TOL)
    {
        return None;
    }

    // 6. coords: √λ · eigenvector for positive eigenvalues; 0 for ~zero; random for negative
    let mut coords = vec![0.0f64; n * dim];
    for a in 0..dim.min(n) {
        let lambda = eig.eigenvalues[order[a]];
        let vcol = eig.eigenvectors.column(order[a]);
        for i in 0..n {
            coords[i * dim + a] = if lambda > EIGVAL_TOL {
                lambda.sqrt() * vcol[i]
            } else if lambda < -EIGVAL_TOL {
                1.0 - 2.0 * rng.gen::<f64>()
            } else {
                0.0
            };
        }
    }
    Some(coords)
}

/// Metric-matrix eigenvalue initialization — RDKit `pickRandomDistMat` + `computeInitialCoords`.
///
/// Pick a random in-bounds distance matrix, form the metric (Gram) matrix
/// `T_ij = ½(d_i0² + d_j0² − d_ij²)` (`d_i0²` = squared distance of atom `i` from the centroid),
/// eigendecompose it, and place atoms with the top-`dim` eigenvectors scaled by `√λ`. An eigenvalue
/// within 0.001 of zero contributes a zero coordinate; a negative one a random coordinate in
/// `[-1, 1]`. Returns `None` on the `computeInitialCoords` reject (see [`coords_from_sq`]) so the
/// caller re-draws, matching RDKit's embed loop; `Some(n_atoms*dim coords)` otherwise.
pub fn metric_matrix(spec: &MoleculeSpec, dim: usize, rng: &mut impl Rng) -> Option<Vec<f64>> {
    let n = spec.n_atoms;

    // 1. random in-bounds distance matrix, squared (symmetric, zero diagonal)
    let mut sq = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..i {
            let (lb, ub) = (spec.lb64(i, j), spec.ub64(i, j));
            let d = lb + rng.gen::<f64>() * (ub - lb);
            let d2 = d * d;
            sq[i * n + j] = d2;
            sq[j * n + i] = d2;
        }
    }

    // 2–6. exact metric embedding of that distance matrix
    coords_from_sq(&sq, n, dim, rng)
}
