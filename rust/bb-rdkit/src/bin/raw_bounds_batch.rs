//! bb-raw-bounds-batch — RDKit's raw (pre-smoothing) setTopolBounds matrix per molecule.
//!
//!   bb-raw-bounds-batch <corpus.smi> > raw_bounds.jsonl
//!
//! The oracle capture for the pure-Rust setTopolBounds port. Each line is one molecule's matrix
//! straight out of `setTopolBounds`, before `triangleSmoothBounds` — the values our Rust constructor
//! must reproduce. One process amortizes RDKit's shared-object load over the whole corpus (a
//! per-process `bb-spec` spends ~1.6 s just linking).
//!
//! Output is JSONL, one object per corpus line:
//!   {"name": "...", "n": <atoms>, "bounds": [<n*n row-major f64, upper-tri=UB lower-tri=LB>]}
//! A molecule RDKit cannot parse is emitted as {"name": "...", "error": "unparseable"} rather than
//! dropped, so the fixture's molecule count matches the corpus.

use std::io::{BufRead, BufWriter, Write};
use std::process::ExitCode;

fn main() -> ExitCode {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: bb-raw-bounds-batch <corpus.smi>");
        return ExitCode::from(2);
    };
    let Ok(f) = std::fs::File::open(&path) else {
        eprintln!("bb-raw-bounds-batch: cannot open {path}");
        return ExitCode::from(2);
    };

    let stdout = std::io::stdout();
    let mut w = BufWriter::new(stdout.lock());
    let (mut ok, mut bad) = (0usize, 0usize);

    for line in std::io::BufReader::new(f).lines() {
        let Ok(line) = line else { break };
        let mut it = line.split_whitespace();
        let (Some(smiles), Some(name)) = (it.next(), it.next()) else {
            continue;
        };
        match bb_rdkit::raw_bounds(smiles) {
            Ok(b) if !b.is_empty() => {
                let n = (b.len() as f64).sqrt() as usize;
                debug_assert_eq!(n * n, b.len(), "raw_bounds returned a non-square matrix");
                // full round-trippable f64 precision: the fixture is the oracle the constructor gate
                // compares against, so it must not itself round at the tolerance scale
                let flat: Vec<String> = b.iter().map(|x| format!("{x}")).collect();
                writeln!(
                    w,
                    r#"{{"name":"{name}","n":{n},"bounds":[{}]}}"#,
                    flat.join(",")
                )
                .expect("write");
                ok += 1;
            }
            // empty vector = RDKit could not parse; Err = RDKit threw
            _ => {
                writeln!(w, r#"{{"name":"{name}","error":"unparseable"}}"#).expect("write");
                bad += 1;
            }
        }
    }
    w.flush().expect("flush");
    eprintln!("  {ok} raw-bounds matrices captured, {bad} unparseable");
    ExitCode::SUCCESS
}
