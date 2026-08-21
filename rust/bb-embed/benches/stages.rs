//! Per-stage benchmarks of the single-conformer embed. `report_counts` prints per-stage wall time for
//! one embed; criterion then times each stage.
//!
//! `BB_SPEC` is required and names the MoleculeSpec JSON to benchmark (`bb-spec-native` produces one):
//!
//!   BB_SPEC=/path/to/spec.json cargo bench -p bb-embed
//!
//! The FF-eval COUNTS in `report_counts` are gated behind the `profile` feature (zero-cost otherwise),
//! so they read 0 unless run as `cargo bench -p bb-embed --features profile`. The wall times are always
//! real. Criterion's "change/regressed" line compares to its own saved baseline (per spec) — meaningless
//! across a spec change, so read the absolute times, not the delta, unless you re-baselined on this spec.

use std::time::Instant;

use bb_core::MoleculeSpec;
use bb_embed::forcefield::{
    self, build_stage_c_constraints, ff_evals, reset_ff_evals, stage_a_energy_grad,
    stage_c_energy_grad, BASIN_DEFAULT,
};
use bb_embed::{init, minimize};
use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use rand::rngs::StdRng;
use rand::SeedableRng;

fn load_spec() -> MoleculeSpec {
    let path = std::env::var("BB_SPEC").unwrap_or_else(|_| {
        panic!(
            "BB_SPEC is not set. Point it at a MoleculeSpec JSON, e.g.\n    \
             bb-spec \"<smiles>\" > /tmp/spec.json\n    \
             BB_SPEC=/tmp/spec.json cargo bench -p bb-embed"
        )
    });
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {path}: {e}"))
}

/// Metric-matrix init, re-drawing on a `computeInitialCoords` reject exactly as the real embed loop
/// does (a `None` is a degenerate sampled distance matrix, not a bad spec — a fresh draw from the same
/// RNG stream usually succeeds). Deterministic given the seed; panics only if 1000 draws all reject.
fn init4(spec: &MoleculeSpec, rng: &mut StdRng) -> Vec<f64> {
    for _ in 0..1000 {
        if let Some(v) = init::metric_matrix(spec, 4, rng) {
            return v;
        }
    }
    panic!("metric_matrix rejected the bench spec on 1000 consecutive draws");
}

fn project_4d_to_3d(a4: &[f64], n: usize) -> Vec<f64> {
    let mut c3 = Vec::with_capacity(n * 3);
    for atom in 0..n {
        for c in 0..3 {
            c3.push(a4[atom * 4 + c]);
        }
    }
    c3
}

/// One seed-path embed, printing FF-eval counts and wall time per stage.
fn report_counts(spec: &MoleculeSpec) {
    let n = spec.n_atoms;
    let mut rng = StdRng::seed_from_u64(1);
    let a0 = init4(spec, &mut rng);

    reset_ff_evals();
    let t = Instant::now();
    let a4 = minimize::minimize_stage_a(spec, a0, 4, &[], BASIN_DEFAULT, 400);
    let (ea, ta) = (ff_evals(), t.elapsed().as_secs_f64() * 1e3);

    reset_ff_evals();
    let t = Instant::now();
    let b4 = minimize::minimize_stage_b(spec, a4, 4, &[], BASIN_DEFAULT, 200);
    let (eb, tb) = (ff_evals(), t.elapsed().as_secs_f64() * 1e3);

    let c3 = project_4d_to_3d(&b4, n);
    let t = Instant::now();
    let (dist_c, angle_c) = build_stage_c_constraints(spec, &c3);
    let tbuild = t.elapsed().as_secs_f64() * 1e3;

    reset_ff_evals();
    let t = Instant::now();
    let final_c = minimize::minimize_stage_c(spec, &dist_c, &angle_c, c3, &[], 300);
    let (ec, tc) = (ff_evals(), t.elapsed().as_secs_f64() * 1e3);

    // Active-set opportunity: how many Stage-C distance constraints are actually violated at the
    // converged geometry? Inactive ones (d within [min,max]) contribute nothing yet are evaluated
    // every iteration — the headroom for a neighbor/active-set optimization.
    let active = dist_c
        .iter()
        .filter(|dc| {
            let d = {
                let (a, b) = (dc.i, dc.j);
                let dx = final_c[3 * a] - final_c[3 * b];
                let dy = final_c[3 * a + 1] - final_c[3 * b + 1];
                let dz = final_c[3 * a + 2] - final_c[3 * b + 2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            d < dc.min_len || d > dc.max_len
        })
        .count();

    // single-eval costs
    let mut rng = StdRng::seed_from_u64(2);
    let x4 = init4(spec, &mut rng);
    reset_ff_evals();
    let t = Instant::now();
    for _ in 0..1000 {
        black_box(stage_a_energy_grad(spec, &x4, 4, BASIN_DEFAULT));
    }
    let ta_one = t.elapsed().as_secs_f64() * 1e3 / 1000.0;

    eprintln!(
        "\n===== per-stage profile ({} atoms, {} dist-constraints, {} angle-constraints) =====",
        n,
        dist_c.len(),
        angle_c.len()
    );
    eprintln!(
        "  Stage A : {ea:5} FF evals, {ta:6.1} ms  ({:.3} ms/eval)",
        ta / ea.max(1) as f64
    );
    eprintln!("  Stage B : {eb:5} FF evals, {tb:6.1} ms",);
    eprintln!("  build C : {:5}            {tbuild:6.1} ms", "-");
    eprintln!(
        "  Stage C : {ec:5} FF evals, {tc:6.1} ms  ({:.3} ms/eval)",
        tc / ec.max(1) as f64
    );
    eprintln!(
        "  total   : {:5} FF evals, {:6.1} ms",
        ea + eb + ec,
        ta + tb + tbuild + tc
    );
    eprintln!("  single stage_a_energy_grad: {ta_one:.4} ms/eval");
    eprintln!(
        "  → Stage-A FF evals ({ea}) vs iter cap (400): {}",
        if ea > 400 {
            "many evals/iter (line search) or 2x cost+grad"
        } else {
            "under cap"
        }
    );
    eprintln!(
        "  ACTIVE-SET: {}/{} dist-constraints active at convergence ({:.0}% inactive → skippable)",
        active,
        dist_c.len(),
        100.0 * (dist_c.len() - active) as f64 / dist_c.len().max(1) as f64
    );
    eprintln!();
}

fn profile(c: &mut Criterion) {
    let spec = load_spec();
    report_counts(&spec);

    let n = spec.n_atoms;
    let mut g = c.benchmark_group("embed_stages");
    g.sample_size(20);

    g.bench_function("metric_matrix_init", |b| {
        let mut rng = StdRng::seed_from_u64(1);
        b.iter(|| black_box(init::metric_matrix(&spec, 4, &mut rng)));
    });
    g.bench_function("stage_a_minimize", |b| {
        let mut rng = StdRng::seed_from_u64(1);
        b.iter_batched(
            || init4(&spec, &mut rng),
            |a| {
                black_box(minimize::minimize_stage_a(
                    &spec,
                    a,
                    4,
                    &[],
                    BASIN_DEFAULT,
                    400,
                ))
            },
            BatchSize::SmallInput,
        );
    });
    g.bench_function("stage_c_minimize", |b| {
        let mut rng = StdRng::seed_from_u64(1);
        // fixed post-Stage-A/B geometry to build constraints from
        let a4 = minimize::minimize_stage_a(
            &spec,
            init4(&spec, &mut rng),
            4,
            &[],
            BASIN_DEFAULT,
            400,
        );
        let b4 = minimize::minimize_stage_b(&spec, a4, 4, &[], BASIN_DEFAULT, 200);
        let c3 = project_4d_to_3d(&b4, n);
        let (dist_c, angle_c) = build_stage_c_constraints(&spec, &c3);
        b.iter_batched(
            || c3.clone(),
            |cc| {
                black_box(minimize::minimize_stage_c(
                    &spec,
                    &dist_c,
                    &angle_c,
                    cc,
                    &[],
                    300,
                ))
            },
            BatchSize::SmallInput,
        );
    });
    g.bench_function("stage_a_energy_grad_once", |b| {
        let mut rng = StdRng::seed_from_u64(1);
        let x = init4(&spec, &mut rng);
        b.iter(|| black_box(stage_a_energy_grad(&spec, &x, 4, BASIN_DEFAULT)));
    });
    g.bench_function("stage_c_energy_grad_once", |b| {
        let mut rng = StdRng::seed_from_u64(1);
        let a4 = minimize::minimize_stage_a(
            &spec,
            init4(&spec, &mut rng),
            4,
            &[],
            BASIN_DEFAULT,
            400,
        );
        let c3 = project_4d_to_3d(&a4, n);
        let (dist_c, angle_c) = build_stage_c_constraints(&spec, &c3);
        b.iter(|| black_box(stage_c_energy_grad(&spec, &dist_c, &angle_c, &c3)));
    });
    let _ = forcefield::reset_ff_evals;
    g.finish();
}

criterion_group!(benches, profile);
criterion_main!(benches);
