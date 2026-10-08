//! bb-batch — read N MoleculeSpec JSONs, embed all of them (rayon over molecules), and report
//! elapsed time and conformers/second on stderr.
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
    // A spec that cannot be embedded in full is named on stderr; the conformers it did embed count.
    let confs: usize = paths
        .par_iter()
        .zip(&specs)
        .map(|(path, s)| match bb_embed::embed_recipe(s, 210185) {
            Ok(blocks) => blocks.iter().map(|b| b.conformers.len()).sum(),
            Err(e) => {
                eprintln!("bb-batch: {path}: {e}");
                match e {
                    bb_embed::RecipeError::Shortfall(shortfall) => shortfall.embedded,
                    _ => 0,
                }
            }
        })
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
