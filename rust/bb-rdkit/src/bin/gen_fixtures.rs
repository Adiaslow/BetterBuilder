//! Fixture generator for the deterministic-component parity goldens.
//!
//! Runs the LIVE RDKit bridge (bb-rdkit) over a corpus and serializes, per case, both the input
//! geometry and RDKit's output — so native gates can validate against these goldens WITHOUT the
//! RDKit toolchain (see `bb-embed/tests/parity_goldens.rs`). The golden is a cache of the live
//! oracle: it is generated only by this bin, never hand-authored, and the input travels with the
//! output so the byte-identical-input coupling that makes a diff mean "formula difference" is
//! preserved. Regenerate on an RDKit version bump; `rust/gates.sh --verify` byte-checks the committed
//! goldens against a fresh live run so they cannot silently rot.
//!
//!   cargo run -p bb-rdkit --bin gen-fixtures -- <corpus.smi> <out_dir>
//!
//! Writes `<out_dir>/VERSION` (the patched-RDKit version, from BB_RDKIT_ROOT) and one `.jsonl` per
//! gate. Deterministic: geometries come from a fixed LCG keyed by molecule index, so regeneration is
//! reproducible.

use std::io::Write;

use bb_embed::forcefield::{BASIN_DEFAULT, W_CHIRAL, W_FOURTH};
use serde::Serialize;

/// Deterministic geometry source (matches the corpus_parity LCG), range ~[-6, 6].
struct Lcg(u64);
impl Lcg {
    fn f(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 12.0
    }
}

/// One Stage-A gradient case: RDKit's gradient at a specific 4D geometry. The geometry is stored so
/// the native gate evaluates its own force field at the identical point.
#[derive(Serialize)]
struct GradRec<'a> {
    smiles: &'a str,
    coords: Vec<f64>,
    grad: Vec<f64>,
}

/// Patched-RDKit version string, from the BB_RDKIT_ROOT directory name (…/rdkit-<version>).
fn rdkit_version() -> String {
    std::env::var("BB_RDKIT_ROOT")
        .ok()
        .and_then(|p| {
            std::path::Path::new(&p)
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.strip_prefix("rdkit-").unwrap_or(s).to_string())
        })
        .unwrap_or_else(|| "unknown".to_string())
}

fn read_smiles(path: &str) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read corpus {path}: {e}"))
        .lines()
        .filter_map(|l| l.split_whitespace().next().map(str::to_string))
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: gen-fixtures <corpus.smi> <out_dir>");
        std::process::exit(2);
    }
    let (corpus, out_dir) = (&args[1], &args[2]);
    std::fs::create_dir_all(out_dir).expect("create out_dir");

    let version = rdkit_version();
    std::fs::write(format!("{out_dir}/VERSION"), format!("{version}\n")).expect("write VERSION");
    eprintln!("patched RDKit version: {version}");

    let smiles = read_smiles(corpus);

    // --- Stage-A gradient goldens ---
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(format!("{out_dir}/stage_a_grad.jsonl")).expect("create stage_a_grad.jsonl"),
    );
    let (mut written, mut spec_fail, mut oracle_fail) = (0usize, 0usize, 0usize);
    for (i, smi) in smiles.iter().enumerate() {
        let spec = match bb_spec::build_native(smi) {
            Ok(s) => s,
            Err(_) => {
                spec_fail += 1;
                continue;
            }
        };
        let n = spec.n_atoms;
        let mut rng = Lcg(0x9E3779B97F4A7C15 ^ (i as u64).wrapping_mul(2654435761));
        let coords: Vec<f64> = (0..n * 4).map(|_| rng.f()).collect();
        let grad = bb_rdkit::stage_a_ff_grad(smi, &coords, W_CHIRAL, W_FOURTH, BASIN_DEFAULT);
        if grad.len() != n * 4 {
            oracle_fail += 1; // RDKit parse/build failure — do not emit a bogus golden
            continue;
        }
        let rec = GradRec { smiles: smi, coords, grad };
        writeln!(out, "{}", serde_json::to_string(&rec).unwrap()).unwrap();
        written += 1;
    }
    out.flush().unwrap();
    eprintln!("stage_a_grad.jsonl: {written} cases ({spec_fail} spec-fail, {oracle_fail} oracle-fail)");

    // --- Stage-B gradient goldens (minimizeFourthDimension: same FF at chiral 0.2 / fourth 1.0) ---
    let mut outb = std::io::BufWriter::new(
        std::fs::File::create(format!("{out_dir}/stage_b_grad.jsonl")).expect("create stage_b_grad.jsonl"),
    );
    let (mut wb, mut bfail) = (0usize, 0usize);
    for (i, smi) in smiles.iter().enumerate() {
        let spec = match bb_spec::build_native(smi) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let n = spec.n_atoms;
        let mut rng = Lcg(0xD1B54A32D192ED03 ^ (i as u64).wrapping_mul(2654435761));
        let coords: Vec<f64> = (0..n * 4).map(|_| rng.f()).collect();
        // Stage B = constructForceField at weights (0.2, 1.0), which the bridge exposes as
        // stage_a_ff_grad called with those weights.
        let grad = bb_rdkit::stage_a_ff_grad(smi, &coords, 0.2, 1.0, BASIN_DEFAULT);
        if grad.len() != n * 4 {
            bfail += 1;
            continue;
        }
        let rec = GradRec { smiles: smi, coords, grad };
        writeln!(outb, "{}", serde_json::to_string(&rec).unwrap()).unwrap();
        wb += 1;
    }
    outb.flush().unwrap();
    eprintln!("stage_b_grad.jsonl: {wb} cases ({bfail} oracle-fail)");

    // --- Stage-C gradient goldens (at a fixed 3D geometry) ---
    // construct3DForceField (torsions incl. the amide-patched force constants + impropers +
    // constraints) at a deterministic 3D geometry. Fixed input, so it is a valid golden AND the
    // host-bridge reference the in-container Stage-C dumper is certified against.
    let mut outc = std::io::BufWriter::new(
        std::fs::File::create(format!("{out_dir}/stage_c_grad.jsonl")).expect("create stage_c_grad.jsonl"),
    );
    let (mut wc, mut cfail) = (0usize, 0usize);
    for (i, smi) in smiles.iter().enumerate() {
        let spec = match bb_spec::build_native(smi) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let n = spec.n_atoms;
        let mut rng = Lcg(0x2545F4914F6CDD1D ^ (i as u64).wrapping_mul(2654435761));
        let coords: Vec<f64> = (0..n * 3).map(|_| rng.f()).collect();
        let grad = bb_rdkit::stage_c_ff_grad(smi, &coords);
        if grad.len() != n * 3 {
            cfail += 1;
            continue;
        }
        let rec = GradRec { smiles: smi, coords, grad };
        writeln!(outc, "{}", serde_json::to_string(&rec).unwrap()).unwrap();
        wc += 1;
    }
    outc.flush().unwrap();
    eprintln!("stage_c_grad.jsonl: {wc} cases ({cfail} oracle-fail)");

    // --- checks goldens (6 pass/reject bits at a fixed 3D geometry) ---
    let mut outk = std::io::BufWriter::new(
        std::fs::File::create(format!("{out_dir}/checks_grad.jsonl")).expect("create checks_grad.jsonl"),
    );
    let (mut wk, mut kfail) = (0usize, 0usize);
    for (i, smi) in smiles.iter().enumerate() {
        let spec = match bb_spec::build_native(smi) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let n = spec.n_atoms;
        let mut rng = Lcg(0x14057B7EF767814F ^ (i as u64).wrapping_mul(2654435761));
        let coords: Vec<f64> = (0..n * 3).map(|_| rng.f()).collect();
        let bits = bb_rdkit::embed_checks(smi, &coords); // Vec<i32>, 6 bits
        if bits.len() != 6 {
            kfail += 1;
            continue;
        }
        let grad: Vec<f64> = bits.iter().map(|&b| b as f64).collect();
        let rec = GradRec { smiles: smi, coords, grad };
        writeln!(outk, "{}", serde_json::to_string(&rec).unwrap()).unwrap();
        wk += 1;
    }
    outk.flush().unwrap();
    eprintln!("checks_grad.jsonl: {wk} cases ({kfail} oracle-fail)");

    // --- bounds goldens (plain smoothed bounds = coord_map with k=0 pins) ---
    let mut outbn = std::io::BufWriter::new(
        std::fs::File::create(format!("{out_dir}/bounds_grad.jsonl")).expect("create bounds_grad.jsonl"),
    );
    let (mut wbn, mut bnfail) = (0usize, 0usize);
    for smi in smiles.iter() {
        let spec = match bb_spec::build_native(smi) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let n = spec.n_atoms;
        let grad = bb_rdkit::coord_map_bounds(smi, &[], &[]); // n*n smoothed bounds, no pins
        if grad.len() != n * n {
            bnfail += 1;
            continue;
        }
        let rec = GradRec { smiles: smi, coords: vec![0.0], grad }; // k=0
        writeln!(outbn, "{}", serde_json::to_string(&rec).unwrap()).unwrap();
        wbn += 1;
    }
    outbn.flush().unwrap();
    eprintln!("bounds_grad.jsonl: {wbn} cases ({bnfail} oracle-fail)");

    // --- coord_map goldens (adjustBoundsMatFromCoordMap): pin the recipe's pin_atoms at an embedded
    //     conformer's positions (valid pins), as corpus_parity does. coords = [k, idx,x,y,z ...]. ---
    let mut outcm = std::io::BufWriter::new(
        std::fs::File::create(format!("{out_dir}/coord_map_grad.jsonl")).expect("create coord_map_grad.jsonl"),
    );
    let (mut wcm, mut cmfail) = (0usize, 0usize);
    for (i, smi) in smiles.iter().enumerate() {
        let spec = match bb_spec::build_native(smi) {
            Ok(s) => s,
            Err(_) => continue,
        };
        if spec.pin_atoms.is_empty() {
            continue;
        }
        let n = spec.n_atoms;
        let confs = match bb_embed::embed(&spec, 1, 0xC0FFEE ^ i as u64) {
            Ok(confs) => confs,
            Err(_) => {
                cmfail += 1;
                continue;
            }
        };
        let conf = &confs[0].coords;
        let (mut idx, mut xyz) = (Vec::new(), Vec::new());
        let mut coords: Vec<f64> = vec![spec.pin_atoms.len() as f64];
        for &a in &spec.pin_atoms {
            let a = a as usize;
            let p = [conf[a * 3], conf[a * 3 + 1], conf[a * 3 + 2]];
            idx.push(a as i32);
            xyz.extend_from_slice(&p);
            coords.push(a as f64);
            coords.extend_from_slice(&p);
        }
        let grad = bb_rdkit::coord_map_bounds(smi, &idx, &xyz);
        if grad.len() != n * n {
            cmfail += 1;
            continue;
        }
        let rec = GradRec { smiles: smi, coords, grad };
        writeln!(outcm, "{}", serde_json::to_string(&rec).unwrap()).unwrap();
        wcm += 1;
    }
    outcm.flush().unwrap();
    eprintln!("coord_map_grad.jsonl: {wcm} cases ({cmfail} embed/oracle-fail)");
}
