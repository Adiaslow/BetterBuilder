//! Corpus-scale parity for the deterministic embed-engine gates — Stage-A/Stage-C force-field
//! gradients, the six post-embed checks, and the core-pin bounds tightening — run native-vs-RDKit
//! over a whole provenanced corpus rather than the hand-picked molecules the per-feature tests use.
//! Those small tests prove each formula on a curated handful; this proves it on the workload, the
//! way `gate_spec` does for the MoleculeSpec (RIGOR.md, Phase 1).
//!
//! Ignored by default (reported as skipped); run with `--ignored` and `BB_PARITY_CORPUS` set to a
//! `<smiles> <name>` file. Every molecule is either compared on every axis that applies to it or
//! named as a failure: both sides rejecting a molecule is agreement, one side alone rejecting it,
//! our embed producing no conformer, or RDKit returning no value are failures, never silent passes.
//! Reports worst-case magnitude per axis and fails if any molecule exceeds tolerance.
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

/// Worst per-component |Δ| relative to gradient magnitude; `None` if RDKit's vector is missing or
/// mis-sized (a bridge failure), which the caller reports as a molecule it could not compare.
fn rel(a: &[f64], b: &[f64]) -> Option<f64> {
    if b.is_empty() || a.len() != b.len() {
        return None;
    }
    let mag = a.iter().chain(b).fold(0.0f64, |m, &x| m.max(x.abs())).max(1.0);
    Some(a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max) / mag)
}

/// One comparison axis over the corpus: molecules compared, the worst difference and where, how many
/// exceeded tolerance, molecules it could not compare (each a failure), and molecules it does not
/// apply to.
struct Axis {
    name: &'static str,
    tol: f64,
    compared: usize,
    worst: f64,
    worst_smi: String,
    over: usize,
    uncompared: Vec<String>,
    not_applicable: usize,
}

impl Axis {
    fn new(name: &'static str, tol: f64) -> Self {
        Axis { name, tol, compared: 0, worst: 0.0, worst_smi: String::new(), over: 0, uncompared: Vec::new(), not_applicable: 0 }
    }
    fn record(&mut self, smi: &str, d: f64) {
        self.compared += 1;
        if d > self.worst {
            self.worst = d;
            self.worst_smi = smi.to_string();
        }
        if d > self.tol {
            self.over += 1;
        }
    }
    fn cannot(&mut self, smi: &str, why: &str) {
        self.uncompared.push(format!("{smi}: {why}"));
    }
    fn report(&self) {
        println!(
            "{:<17}: {} compared, worst {:.2e} ({}) | {} over {:.0e} | {} could not compare | {} not applicable",
            self.name, self.compared, self.worst, self.worst_smi, self.over, self.tol, self.uncompared.len(), self.not_applicable
        );
        for u in self.uncompared.iter().take(8) {
            println!("    {u}");
        }
    }
    fn failures(&self) -> usize {
        self.over + self.uncompared.len()
    }
}

// Fixed tolerances, not yet derived from the evaluations' floating-point error.
const TOL_A: f64 = 1e-12; // f64 bounds → Stage-A distance term matches RDKit to machine precision
const TOL_B: f64 = 1e-12; // Stage B = the same distance term at weights 0.2/1.0
const TOL_C: f64 = 1e-6;  // machine-precision except the FP floor at V=100 planarity torsions (near-planar 1/sinφ)
const TOL_CMB: f64 = 1e-9; // f64 raw_bounds + f64 coordMap → machine precision
const TOL_RE: f64 = 1e-12; // reject energy is an exact identity (native + chiral+4th) vs RDKit calcEnergy; only summation order differs

#[test]
#[ignore = "corpus-scale gate: needs BB_PARITY_CORPUS; run with --ignored, as rust/gates.sh live and hunt do"]
fn corpus_parity() {
    let corpus = std::env::var("BB_PARITY_CORPUS").expect("BB_PARITY_CORPUS must name the corpus to check");
    let text = std::fs::read_to_string(&corpus).expect("read corpus");
    let smiles: Vec<String> = text
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_string))
        .collect();

    let mut a = Axis::new("Stage-A gradient", TOL_A);
    let mut b = Axis::new("Stage-B gradient", TOL_B);
    let mut c = Axis::new("Stage-C gradient", TOL_C);
    let mut re = Axis::new("reject energy", TOL_RE);
    let mut cmb = Axis::new("coord_map_bounds", TOL_CMB);
    // the checks' accept/reject decision: a disagreement counts as 1.0 against a tolerance of 0
    let mut dec = Axis::new("checks decision", 0.0);
    let (mut n, mut both_reject) = (0usize, 0usize);
    let (mut only_ours_rejects, mut only_rdkit_rejects, mut no_conformer) = (Vec::new(), Vec::new(), Vec::new());

    for (i, smi) in smiles.iter().enumerate() {
        // setup: both sides must agree on whether the molecule is acceptable at all
        let ours = bb_spec::build_native(smi);
        let theirs = bb_rdkit::build_spec(smi);
        let spec = match (ours, theirs) {
            (Ok(s), Ok(_)) => s,
            (Err(_), Err(_)) => {
                both_reject += 1;
                continue;
            }
            (Err(e), Ok(_)) => {
                only_ours_rejects.push(format!("{smi}: {e}"));
                continue;
            }
            (Ok(_), Err(e)) => {
                only_rdkit_rejects.push(format!("{smi}: {e}"));
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
        match rel(&na, &ra) {
            Some(d) => a.record(smi, d),
            None => a.cannot(smi, "RDKit returned no Stage-A gradient"),
        }

        // --- Stage-A reject ENERGY at the same geometry: native's `stage_a_reject_energy` must equal
        // RDKit's `calcEnergy` (the quantity firstMinimization thresholds at 0.05/atom). This random
        // geometry activates all three terms (distance + chiral-volume + 4th-dim), so it genuinely
        // exercises the RDKit-convention (w·x², no ½) chiral/4th energy the reject depends on. If the
        // energy matches to machine precision the 0.05 reject decision is identical by construction. ---
        let nre = forcefield::stage_a_reject_energy(&spec, &coords4, 4, forcefield::BASIN_DEFAULT);
        let rre = bb_rdkit::stage_a_energy(smi, &coords4, forcefield::W_CHIRAL, forcefield::W_FOURTH, forcefield::BASIN_DEFAULT);
        if rre.is_nan() {
            re.cannot(smi, "RDKit returned NaN for the reject energy");
        } else {
            re.record(smi, (nre - rre).abs() / nre.abs().max(rre.abs()).max(1.0));
        }

        // --- Stage B: same distance term at chiral 0.2 / 4th-dim 1.0 (minimizeFourthDimension) ---
        let nb2 = forcefield::stage_b_energy_grad(&spec, &coords4, 4, forcefield::BASIN_DEFAULT).1;
        let rb2 = bb_rdkit::stage_a_ff_grad(smi, &coords4, 0.2, 1.0, forcefield::BASIN_DEFAULT);
        match rel(&nb2, &rb2) {
            Some(d) => b.record(smi, d),
            None => b.cannot(smi, "RDKit returned no Stage-B gradient"),
        }

        // --- Everything geometry-based needs one embedded conformer ---
        let confs = match bb_embed::embed(&spec, 1, 0xC0FFEE ^ i as u64) {
            Ok(confs) => confs,
            Err(e) => {
                no_conformer.push(format!("{smi}: {e}"));
                continue;
            }
        };
        let conf = &confs[0].coords;

        // --- Stage C full gradient ---
        let (dc, ac) = forcefield::build_stage_c_constraints(&spec, conf);
        let nc = forcefield::stage_c_energy_grad(&spec, &dc, &ac, conf).1;
        let rc = bb_rdkit::stage_c_ff_grad(smi, conf);
        match rel(&nc, &rc) {
            Some(d) => c.record(smi, d),
            None => c.cannot(smi, "RDKit returned no Stage-C gradient"),
        }

        // --- The checks' accept/reject decision, native vs RDKit's real routines, on the same
        // conformer: RDKit accepts only if every check passes (`embed_checks_parity` covers each
        // check deciding alone) ---
        let rbits = bb_rdkit::embed_checks(smi, conf);
        if rbits.len() == 6 {
            let rdkit_accepts = rbits.iter().all(|&x| x == 1);
            let native_accepts = checks::passes_checks(&spec, conf, &[]);
            dec.record(smi, if native_accepts == rdkit_accepts { 0.0 } else { 1.0 });
        } else {
            dec.cannot(smi, "RDKit returned no check results");
        }

        // --- Core-pin bounds tightening with the recipe's real pins ---
        if spec.pin_atoms.is_empty() {
            cmb.not_applicable += 1;
        } else {
            let mut pinned = vec![None; np];
            let (mut idx, mut xyz) = (Vec::new(), Vec::new());
            for &p in &spec.pin_atoms {
                let p = p as usize;
                let x = [conf[p * 3], conf[p * 3 + 1], conf[p * 3 + 2]];
                pinned[p] = Some(x);
                idx.push(p as i32);
                xyz.extend_from_slice(&x);
            }
            let ncmb = bb_embed::bounds::coord_map_bounds_f64(&spec.raw_bounds_f64, np, &pinned);
            let rcmb = bb_rdkit::coord_map_bounds(smi, &idx, &xyz);
            if rcmb.is_empty() || ncmb.len() != rcmb.len() {
                cmb.cannot(smi, "RDKit returned no coordMap bounds of the right size");
            } else {
                cmb.record(smi, ncmb.iter().zip(&rcmb).map(|(&x, &y)| (x - y).abs()).fold(0.0f64, f64::max));
            }
        }
    }

    println!("\n=== corpus parity: {corpus} ===");
    println!(
        "molecules: {} in corpus, {n} compared, {both_reject} rejected by both sides, {} rejected by ours only, {} by RDKit only, {} without a conformer",
        smiles.len(), only_ours_rejects.len(), only_rdkit_rejects.len(), no_conformer.len()
    );
    for e in only_ours_rejects.iter().chain(&only_rdkit_rejects).chain(&no_conformer).take(8) {
        println!("    {e}");
    }
    for axis in [&a, &b, &c, &re, &dec, &cmb] {
        axis.report();
    }

    assert!(n > 0, "no molecule of {} was compared", smiles.len());
    assert!(only_ours_rejects.is_empty(), "{} molecules rejected by our setup only", only_ours_rejects.len());
    assert!(only_rdkit_rejects.is_empty(), "{} molecules rejected by RDKit only", only_rdkit_rejects.len());
    assert!(no_conformer.is_empty(), "{} molecules got no conformer to compare on", no_conformer.len());
    for axis in [&a, &b, &c, &re, &dec, &cmb] {
        assert_eq!(axis.failures(), 0, "{}: {} over {:.0e}, {} could not compare", axis.name, axis.over, axis.tol, axis.uncompared.len());
    }
}
