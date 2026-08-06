//! Initial coordinates for the DistGeom embed. Two modes, matching RDKit:
//! - [`random_box`] — `useRandomCoords=True` (the sidechain embeds).
//! - [`metric_matrix`] — `useRandomCoords=False` (the seed/core embeds): a distance-geometry-informed
//!   metric-matrix eigenvalue embedding, so seeds cluster instead of scattering.

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

/// Metric-matrix eigenvalue initialization — RDKit `pickRandomDistMat` + `computeInitialCoords`.
///
/// Pick a random in-bounds distance matrix, form the metric (Gram) matrix
/// `T_ij = ½(d_i0² + d_j0² − d_ij²)` (`d_i0²` = squared distance of atom `i` from the centroid),
/// eigendecompose it, and place atoms with the top-`dim` eigenvectors scaled by `√λ`. Because the
/// distances respect the bounds, the resulting structure is already near-valid, so conformers
/// cluster (unlike the random box). Returns `n_atoms*dim` coords.
pub fn metric_matrix(spec: &MoleculeSpec, dim: usize, rng: &mut impl Rng) -> Vec<f64> {
    let n = spec.n_atoms;

    // 1. random in-bounds distance matrix, squared (symmetric, zero diagonal)
    let mut sq = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..i {
            let (lb, ub) = (spec.lb(i, j) as f64, spec.ub(i, j) as f64);
            let d = lb + rng.gen::<f64>() * (ub - lb);
            let d2 = d * d;
            sq[i * n + j] = d2;
            sq[j * n + i] = d2;
        }
    }

    // 2. sumSqD2 = (Σ_{i>j} d_ij²) / N²  (lower triangle only, matching RDKit's SymmMatrix data)
    let mut sum_lower = 0.0;
    for i in 0..n {
        for j in 0..i {
            sum_lower += sq[i * n + j];
        }
    }
    let sum_sq_d2 = sum_lower / (n * n) as f64;

    // 3. squared distances from the centroid: d_i0² = (Σ_j d_ij²)/N − sumSqD2
    let mut d0 = vec![0.0f64; n];
    for i in 0..n {
        let row: f64 = (0..n).map(|j| sq[i * n + j]).sum();
        d0[i] = row / n as f64 - sum_sq_d2;
    }

    // 4. metric (Gram) matrix T_ij = ½(d_i0² + d_j0² − d_ij²)
    let mut t = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            t[(i, j)] = 0.5 * (d0[i] + d0[j] - sq[i * n + j]);
        }
    }

    // 5. symmetric eigendecomposition; take the top `dim` eigenpairs (largest eigenvalue first)
    let eig = t.symmetric_eigen();
    let mut order: Vec<usize> = (0..n).collect();
    // Descending eigenvalue order. `total_cmp` gives a total order with no panic path (the eigenvalues
    // of a symmetric matrix are real, so this never sees a NaN, but total_cmp needs no such assumption).
    order.sort_by(|&a, &b| eig.eigenvalues[b].total_cmp(&eig.eigenvalues[a]));

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
    coords
}
