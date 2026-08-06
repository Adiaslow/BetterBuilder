//! Batch throughput benchmark: read N MoleculeSpec JSONs, embed all (rayon over molecules), report
//! conformers/second. This is the production pattern — the workload we'd fan out across cores.
//!
//!   bb-batch <spec1.json> <spec2.json> ...        (RAYON_NUM_THREADS controls core count)

use std::time::Instant;

use bb_core::MoleculeSpec;
use rayon::prelude::*;

fn main() {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        eprintln!("usage: bb-batch <spec.json> ...");
        std::process::exit(1);
    }
    let specs: Vec<MoleculeSpec> = paths
        .iter()
        .map(|p| {
            serde_json::from_str(&std::fs::read_to_string(p).expect("read spec"))
                .expect("parse spec")
        })
        .collect();

    let threads = rayon::current_num_threads();
    let t = Instant::now();
    // Molecule-level parallelism; embed_recipe also parallelizes internally (one shared rayon pool).
    let confs: usize = specs
        .par_iter()
        .map(|s| bb_embed::embed_recipe(s, 210185).len())
        .sum();
    let dt = t.elapsed().as_secs_f64();

    eprintln!(
        "{} molecules, {} conformers, {} threads: {:.2}s  →  {:.0} conf/s  ({:.1} ms/conf-core)",
        specs.len(),
        confs,
        threads,
        dt,
        confs as f64 / dt,
        dt * threads as f64 / confs as f64 * 1e3,
    );
}
