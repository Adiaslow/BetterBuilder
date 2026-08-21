//! Correctness of the metric-matrix initial embedding ([`bb_embed::init::coords_from_sq`], and the
//! [`metric_gram`](bb_embed::init::metric_gram) Gram-matrix construction it uses).
//!
//! Ground truth is independent of the implementation: take a random point cloud in `d` dimensions,
//! form its EXACT pairwise squared-distance matrix, embed that matrix with native's metric-matrix
//! method in the same `d`, and require the embedding to reproduce the original inter-point distances
//! to machine precision. A `d`-embeddable distance matrix has a rank-`d` Gram matrix (its top `d`
//! eigenvalues positive, the rest ~0), so classical MDS recovers the geometry exactly up to a rigid
//! motion — which the distance matrix is invariant to. Any error in the double-centering formula,
//! the eigenvalue selection, or the `√λ·v` placement breaks the reconstruction and fails the test.
//!
//! This is a native-correctness gate (no RDKit oracle): native deliberately uses an exact
//! eigensolver, not RDKit's power iteration, so there is no parity target here — the source of truth
//! is the geometry itself.

use bb_embed::init;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Exact pairwise squared-distance matrix (`n*n`) of `n` points laid out `n*d` row-major.
fn sq_dist(points: &[f64], n: usize, d: usize) -> Vec<f64> {
    let mut sq = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for c in 0..d {
                let diff = points[i * d + c] - points[j * d + c];
                s += diff * diff;
            }
            sq[i * n + j] = s;
        }
    }
    sq
}

/// Recover distances (not squared) from an embedding and compare to the target, worst |Δ| in Å.
fn worst_distance_error(embed: &[f64], target_sq: &[f64], n: usize, d: usize) -> f64 {
    let got = sq_dist(embed, n, d);
    let mut worst = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let a = got[i * n + j].max(0.0).sqrt();
            let b = target_sq[i * n + j].max(0.0).sqrt();
            worst = worst.max((a - b).abs());
        }
    }
    worst
}

#[test]
fn recovers_embeddable_geometry_to_machine_precision() {
    // A d-embeddable matrix must be reconstructed exactly (up to rigid motion). Threshold is a
    // machine-precision-scale bound: the reconstruction is exact arithmetic on an exactly rank-d
    // Gram matrix, so residual is float round-off, not a fitted tolerance.
    const TOL: f64 = 1e-9;

    // deterministic seed (fixed constant; no entropy/time dependence)
    let mut rng = StdRng::seed_from_u64(0x9E37_79B9_7F4A_7C15);

    let mut worst_overall = 0.0f64;
    let mut cases = 0usize;
    for &d in &[3usize, 4] {
        for &n in &[4usize, 5, 7, 10, 16, 24] {
            if n <= d {
                continue;
            }
            for trial in 0..8 {
                // random point cloud in d dimensions, spread over ~[-5,5]
                let points: Vec<f64> = (0..n * d).map(|_| 10.0 * (rng.gen::<f64>() - 0.5)).collect();
                let target_sq = sq_dist(&points, n, d);

                // A d-embeddable matrix is non-degenerate (top-d eigenvalues positive, sqD0i > 0), so
                // the computeInitialCoords reject never fires here — coords_from_sq returns Some.
                let embed = init::coords_from_sq(&target_sq, n, d, &mut rng)
                    .expect("embeddable geometry must not trigger the init reject");
                assert_eq!(embed.len(), n * d, "embedding wrong length");
                let worst = worst_distance_error(&embed, &target_sq, n, d);
                worst_overall = worst_overall.max(worst);
                cases += 1;
                assert!(
                    worst < TOL,
                    "d={d} n={n} trial={trial}: distance reconstruction error {worst:.3e} Å exceeds {TOL:.0e}"
                );
            }
        }
    }
    // Non-vacuity: the loop must actually have exercised reconstructions.
    assert!(cases >= 40, "too few cases exercised: {cases}");
    println!("embedding recovery: {cases} embeddable geometries, worst distance error {worst_overall:.3e} Å");
}
