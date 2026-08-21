//! Deterministic-component parity, validated against committed RDKit goldens — **no RDKit toolchain
//! required**. The goldens under `validation/fixtures/rdkit/` are a cache of the live bb-rdkit bridge
//! oracle (generated only by `cargo run -p bb-rdkit --bin gen-fixtures`); each record carries the
//! input geometry alongside RDKit's output, so native evaluates its own force field at the *identical*
//! point and a difference means a formula difference, exactly as in the live comparison.
//!
//! This is the fast inner-loop gate. The authoritative anti-staleness check is `gates.sh --verify`,
//! which regenerates from live RDKit and byte-diffs these committed goldens (pinned to VERSION).
//!
//! Missing goldens are a hard failure (not a skip): a parity gate with no reference is vacuous.

use bb_embed::forcefield::{stage_a_energy_grad, BASIN_DEFAULT};
use serde::Deserialize;

/// Repo-root-anchored fixture path (independent of the test binary's working directory).
fn fixture(rel: &str) -> String {
    format!("{}/../../validation/fixtures/{}", env!("CARGO_MANIFEST_DIR"), rel)
}

#[derive(Deserialize)]
struct GradRec {
    smiles: String,
    coords: Vec<f64>,
    grad: Vec<f64>,
}

/// Worst per-component |Δ| relative to gradient magnitude (matches ff_parity's `rel_grad_diff`).
fn rel_grad(a: &[f64], b: &[f64]) -> f64 {
    let mag = a.iter().chain(b).fold(0.0f64, |m, &x| m.max(x.abs())).max(1.0);
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max) / mag
}

#[test]
fn stage_a_gradient_matches_golden() {
    // f64 bounds → Stage-A distance term matches RDKit to machine precision (RIGOR.md).
    const TOL: f64 = 1e-12;

    let path = fixture("rdkit/stage_a_grad.jsonl");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!("golden missing ({path}): {e}\n  regenerate: cargo run -p bb-rdkit --bin gen-fixtures -- validation/seeds_100.smi validation/fixtures/rdkit")
    });

    let (mut n, mut fails, mut worst) = (0usize, 0usize, 0.0f64);
    let mut worst_smi = String::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: GradRec = serde_json::from_str(line).expect("parse golden record");
        let spec = bb_spec::build_native(&rec.smiles).expect("native spec");
        assert_eq!(
            rec.coords.len(),
            spec.n_atoms * 4,
            "{}: golden geometry length {} != native n_atoms*4 {}",
            rec.smiles,
            rec.coords.len(),
            spec.n_atoms * 4
        );
        let native = stage_a_energy_grad(&spec, &rec.coords, 4, BASIN_DEFAULT).1;
        let r = rel_grad(&native, &rec.grad);
        if r > worst {
            worst = r;
            worst_smi = rec.smiles.clone();
        }
        if r > TOL {
            fails += 1;
        }
        n += 1;
    }

    println!("stage_a_gradient vs golden: {n} cases, worst rel {worst:.2e} ({worst_smi}), {fails} over {TOL:.0e}");
    assert!(n > 0, "no golden cases — empty fixture");
    assert_eq!(fails, 0, "Stage-A gradient exceeded {TOL:.0e} vs golden on {fails} molecules");
}
