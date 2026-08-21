//! Time torsion assignment with the same discipline as the C++ baseline.
//!
//!   bench-torsions <corpus.smi>
//!
//! `validation/bench/time_cpp_torsions.cc` times C++ RDKit doing SMILES -> AddHs ->
//! getExperimentalTorsions over a corpus, after a warm-up and excluding process startup and I/O.
//! This mirrors that exactly for the pure-Rust path — read the corpus up front, compile the library
//! once, warm up on the first few molecules, then time only the per-molecule perception plus
//! assignment. Same work, same exclusions, so the two ms/mol numbers compare like for like.
//!
//! It reports the total torsion count too: it must equal the C++ harness's count, or the two are
//! not doing the same work and the timing is meaningless.

use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("usage: bench-torsions <corpus.smi>");
    let text = std::fs::read_to_string(&path).expect("read corpus");
    let smis: Vec<&str> = text
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .collect();

    let lib = bb_perceive::torsion_lib::patterns();
    let asts = bb_perceive::torsions::compile(&lib);

    // returns (parsed ok, torsion count) in one pass, so the timed loop perceives each molecule once
    let run_one = |smi: &str| -> (bool, usize) {
        match bb_perceive::smarts_match::perceive(smi) {
            Ok(mol) => {
                let (mut ts, mut done) = bb_perceive::torsions::assign(&mol, &lib, &asts);
                let (bk, _imp) = bb_perceive::torsions::basic_knowledge(&mol, &mut done);
                ts.extend(bk);
                (true, ts.len())
            }
            Err(_) => (false, 0),
        }
    };

    // warm-up, matching the C++ harness's first-five warm-up
    for smi in smis.iter().take(5) {
        run_one(smi);
    }

    let mut n_ok = 0usize;
    let mut n_tors = 0usize;
    let t0 = Instant::now();
    for smi in &smis {
        let (ok, c) = run_one(smi);
        if ok {
            n_ok += 1;
        }
        n_tors += c;
    }
    let secs = t0.elapsed().as_secs_f64();

    eprintln!(
        "rust bb-perceive e2e (SMILES->perceive->assign): {n_ok} mols in {secs:.5} s = {:.5} ms/mol, {n_tors} torsions",
        1000.0 * secs / n_ok as f64
    );
}
