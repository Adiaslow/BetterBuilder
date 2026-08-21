//! bb-spec-native-batch — a whole corpus through the RDKit-free assembler in one process.
//!
//!   bb-spec-native-batch <corpus.smi> > specs.jsonl
//!
//! One JSON object per line, keyed by the corpus name — the pure-Rust counterpart to `bb-spec-batch`
//! (the RDKit bridge), so a corpus can be gated JSONL-to-JSONL in two process runs instead of one
//! subprocess launch per molecule. A molecule that cannot be built is emitted as
//! `{"name":..,"error":..}` rather than dropped, so the line count always matches the corpus.

use std::io::{BufRead, BufWriter, Write};
use std::process::ExitCode;

fn main() -> ExitCode {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: bb-spec-native-batch <corpus.smi>");
        return ExitCode::from(2);
    };
    let Ok(f) = std::fs::File::open(&path) else {
        eprintln!("bb-spec-native-batch: cannot open {path}");
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
        match bb_spec::build_native(smiles) {
            Ok(spec) => {
                let body = serde_json::to_string(&spec).expect("serialize spec");
                writeln!(w, r#"{{"name":"{name}","spec":{body}}}"#).expect("write");
                ok += 1;
            }
            Err(e) => {
                let msg = e.replace('"', "'");
                writeln!(w, r#"{{"name":"{name}","error":"{msg}"}}"#).expect("write");
                bad += 1;
            }
        }
    }
    w.flush().expect("flush");
    eprintln!("  {ok} specs built, {bad} rejected");
    ExitCode::SUCCESS
}
