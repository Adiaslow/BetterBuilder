//! Corpus-scale parity for the deterministic embed-engine gates — Stage-A/Stage-C force-field
//! gradients, the six post-embed checks, and the core-pin bounds tightening — run native-vs-RDKit
//! over a whole provenanced corpus rather than the hand-picked molecules the per-feature tests use.
//! Those small tests prove each formula on a curated handful; this proves it on the workload, the
//! way `gate_spec` does for the MoleculeSpec (RIGOR.md, Phase 1).
//!
//! Opt-in: set `BB_PARITY_CORPUS` to a `<smiles> <name>` file. Without it the test skips (so it does
//! not slow the default `cargo test`). Reports worst-case magnitude and a disagreement count per
//! axis, and fails if any molecule exceeds tolerance — the count *is* the result.
//!
//!   BB_PARITY_CORPUS=validation/seeds_100.smi cargo test -p bb-rdkit --test corpus_parity -- --nocapture

use bb_embed::{checks, forcefield};

struct Lcg(u64);
impl Lcg {
    fn f(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 12.0
    }
}

/// Worst per-component |Δ| relative to gradient magnitude; None if the oracle vector is missing or
/// mis-sized (a bridge parse failure), so the caller can count it as an oracle skip, not a pass.
fn rel(a: &[f64], b: &[f64]) -> Option<f64> {
    if b.is_empty() || a.len() != b.len() {
        return None;
    }
    let mag = a.iter().chain(b).fold(0.0f64, |m, &x| m.max(x.abs())).max(1.0);
    Some(a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max) / mag)
}

fn native_check_bits(spec: &bb_core::MoleculeSpec, c: &[f64]) -> [i32; 6] {
    let pinned = vec![None; spec.n_atoms];
    [
        checks::check_tetrahedral_centers(spec, c, &pinned) as i32,
        checks::check_chiral_centers(spec, c) as i32,
        checks::planarity_ok(spec, c) as i32,
        checks::double_bond_geometry_ok(spec, c) as i32,
        checks::final_chiral_checks(spec, c) as i32,
        checks::double_bond_stereo_ok(spec, c) as i32,
    ]
}

// Derived tolerances (RIGOR.md): Stage A ~machine precision; Stage C carries f32-bounds + torsion
// numerics; coord_map_bounds is f32 storage of a f64 matrix. Checks are exact (bit agreement).
const TOL_A: f64 = 1e-12; // f64 bounds → Stage-A distance term matches RDKit to machine precision
const TOL_B: f64 = 1e-12; // Stage B = the same distance term at weights 0.2/1.0
const TOL_C: f64 = 1e-6;  // machine-precision except the FP floor at V=100 planarity torsions (near-planar 1/sinφ)
const TOL_CMB: f64 = 1e-9; // f64 raw_bounds + f64 coordMap → machine precision
const TOL_RE: f64 = 1e-12; // reject energy is an exact identity (native + chiral+4th) vs RDKit calcEnergy; only summation order differs

#[test]
fn corpus_parity() {
    let corpus = match std::env::var("BB_PARITY_CORPUS") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("BB_PARITY_CORPUS unset — skipping corpus-scale parity");
            return;
        }
    };
    let text = std::fs::read_to_string(&corpus).expect("read corpus");
    let smiles: Vec<String> = text
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_string))
        .collect();

    let (mut n, mut spec_fail, mut embed_fail) = (0usize, 0usize, 0usize);
    let (mut worst_a, mut worst_b, mut worst_c, mut worst_cmb) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut fail_a, mut fail_b, mut fail_c, mut fail_cmb) = (0usize, 0usize, 0usize, 0usize);
    let (mut worst_re, mut fail_re) = (0.0f64, 0usize);
    let (mut check_diffs, mut check_evals) = (0usize, 0usize);
    let (mut wa_smi, mut wb_smi, mut wc_smi, mut wcmb_smi) = (String::new(), String::new(), String::new(), String::new());
    let mut wre_smi = String::new();
    let mut check_diff_ex: Vec<String> = Vec::new();

    for (i, smi) in smiles.iter().enumerate() {
        let spec = match bb_spec::build_native(smi) {
            Ok(s) => s,
            Err(_) => {
                spec_fail += 1;
                continue;
            }
        };
        n += 1;
        let np = spec.n_atoms;

        // --- Stage A at a random 4D geometry (formula check, no distance-constraint coupling) ---
        let mut rng = Lcg(0x9E3779B97F4A7C15 ^ (i as u64).wrapping_mul(2654435761));
        let coords4: Vec<f64> = (0..np * 4).map(|_| rng.f()).collect();
        let na = forcefield::stage_a_energy_grad(&spec, &coords4, 4, forcefield::BASIN_DEFAULT).1;
        let ra = bb_rdkit::stage_a_ff_grad(smi, &coords4, forcefield::W_CHIRAL, forcefield::W_FOURTH, forcefield::BASIN_DEFAULT);
        if let Some(r) = rel(&na, &ra) {
            if r > worst_a {
                worst_a = r;
                wa_smi = smi.clone();
            }
            if r > TOL_A {
                fail_a += 1;
            }
        }

        // --- Stage-A reject ENERGY at the same geometry: native's `stage_a_reject_energy` must equal
        // RDKit's `calcEnergy` (the quantity firstMinimization thresholds at 0.05/atom). This random
        // geometry activates all three terms (distance + chiral-volume + 4th-dim), so it genuinely
        // exercises the RDKit-convention (w·x², no ½) chiral/4th energy the reject depends on. If the
        // energy matches to machine precision the 0.05 reject decision is identical by construction. ---
        let nre = forcefield::stage_a_reject_energy(&spec, &coords4, 4, forcefield::BASIN_DEFAULT);
        let rre = bb_rdkit::stage_a_energy(smi, &coords4, forcefield::W_CHIRAL, forcefield::W_FOURTH, forcefield::BASIN_DEFAULT);
        if !rre.is_nan() {
            let r = (nre - rre).abs() / nre.abs().max(rre.abs()).max(1.0);
            if r > worst_re {
                worst_re = r;
                wre_smi = smi.clone();
            }
            if r > TOL_RE {
                fail_re += 1;
            }
        }

        // --- Stage B: same distance term at chiral 0.2 / 4th-dim 1.0 (minimizeFourthDimension) ---
        let nb2 = forcefield::stage_b_energy_grad(&spec, &coords4, 4, forcefield::BASIN_DEFAULT).1;
        let rb2 = bb_rdkit::stage_a_ff_grad(smi, &coords4, 0.2, 1.0, forcefield::BASIN_DEFAULT);
        if let Some(r) = rel(&nb2, &rb2) {
            if r > worst_b {
                worst_b = r;
                wb_smi = smi.clone();
            }
            if r > TOL_B {
                fail_b += 1;
            }
        }

        // --- Everything geometry-based needs one embedded conformer ---
        let confs = bb_embed::embed(&spec, 1, 0xC0FFEE ^ i as u64);
        let conf = match confs.first() {
            Some(c) if c.coords.len() == np * 3 => &c.coords,
            _ => {
                embed_fail += 1;
                continue;
            }
        };

        // --- Stage C full gradient ---
        let (dc, ac) = forcefield::build_stage_c_constraints(&spec, conf);
        let nc = forcefield::stage_c_energy_grad(&spec, &dc, &ac, conf).1;
        let rc = bb_rdkit::stage_c_ff_grad(smi, conf);
        if let Some(r) = rel(&nc, &rc) {
            if r > worst_c {
                worst_c = r;
                wc_smi = smi.clone();
            }
            if r > TOL_C {
                fail_c += 1;
            }
        }

        // --- The six checks, native vs RDKit's real routines, on the same conformer ---
        let rb = bb_rdkit::embed_checks(smi, conf);
        if rb.len() == 6 {
            check_evals += 1;
            let nb = native_check_bits(&spec, conf);
            if nb[..] != rb[..] {
                check_diffs += 1;
                if check_diff_ex.len() < 8 {
                    check_diff_ex.push(format!("{smi}: native={nb:?} rdkit={rb:?}"));
                }
            }
        }

        // --- Core-pin bounds tightening with the recipe's real pins ---
        if !spec.pin_atoms.is_empty() && !spec.raw_bounds_f64.is_empty() {
            let mut pinned = vec![None; np];
            let (mut idx, mut xyz) = (Vec::new(), Vec::new());
            for &a in &spec.pin_atoms {
                let a = a as usize;
                let p = [conf[a * 3], conf[a * 3 + 1], conf[a * 3 + 2]];
                pinned[a] = Some(p);
                idx.push(a as i32);
                xyz.extend_from_slice(&p);
            }
            let ncmb = bb_embed::bounds::coord_map_bounds_f64(&spec.raw_bounds_f64, np, &pinned);
            let rcmb = bb_rdkit::coord_map_bounds(smi, &idx, &xyz);
            if !rcmb.is_empty() && ncmb.len() == rcmb.len() {
                let d = ncmb.iter().zip(&rcmb).map(|(&x, &y)| (x - y).abs()).fold(0.0f64, f64::max);
                if d > worst_cmb {
                    worst_cmb = d;
                    wcmb_smi = smi.clone();
                }
                if d > TOL_CMB {
                    fail_cmb += 1;
                }
            }
        }
    }

    println!("\n=== corpus parity: {corpus} ===");
    println!("molecules: {n} built, {spec_fail} spec-fail, {embed_fail} embed-fail");
    println!("Stage-A gradient : worst {worst_a:.2e} ({wa_smi})  | {fail_a} over {TOL_A:.0e}");
    println!("Stage-B gradient : worst {worst_b:.2e} ({wb_smi})  | {fail_b} over {TOL_B:.0e}");
    println!("Stage-C gradient : worst {worst_c:.2e} ({wc_smi})  | {fail_c} over {TOL_C:.0e}");
    println!("checks (6-bit)   : {check_diffs} disagree / {check_evals} evaluated");
    for e in &check_diff_ex {
        println!("    {e}");
    }
    println!("coord_map_bounds : worst {worst_cmb:.2e} ({wcmb_smi})  | {fail_cmb} over {TOL_CMB:.0e}");
    println!("reject energy    : worst {worst_re:.2e} ({wre_smi})  | {fail_re} over {TOL_RE:.0e}");

    assert_eq!(fail_a, 0, "Stage-A gradient exceeded {TOL_A:.0e} on {fail_a} molecules");
    assert_eq!(fail_b, 0, "Stage-B gradient exceeded {TOL_B:.0e} on {fail_b} molecules");
    assert_eq!(fail_c, 0, "Stage-C gradient exceeded {TOL_C:.0e} on {fail_c} molecules");
    assert_eq!(check_diffs, 0, "checks disagreed on {check_diffs} molecules");
    assert_eq!(fail_cmb, 0, "coord_map_bounds exceeded {TOL_CMB:.0e} on {fail_cmb} molecules");
    assert_eq!(fail_re, 0, "reject energy exceeded {TOL_RE:.0e} on {fail_re} molecules");
}
