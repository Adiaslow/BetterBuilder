//! bb-spec-batch — a whole corpus through the patched RDKit in one process.
//!
//!   bb-spec-batch <corpus.smi> > specs.jsonl
//!
//! One JSON object per line, keyed by the corpus name. `bb-spec` builds one spec per process, which
//! spends its time loading RDKit's shared objects rather than on the molecule; capturing a fixture
//! over thousands of molecules pays that cost once here instead.
//!
//! A molecule the patched RDKit rejects is emitted as `{"name":..,"error":..}` rather than dropped,
//! so a fixture's molecule count always matches the corpus.

use std::io::{BufRead, BufWriter, Write};
use std::process::ExitCode;

fn main() -> ExitCode {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: bb-spec-batch <corpus.smi>");
        return ExitCode::from(2);
    };
    let Ok(f) = std::fs::File::open(&path) else {
        eprintln!("bb-spec-batch: cannot open {path}");
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
        match bb_rdkit::build_spec(smiles) {
            Ok(spec) => {
                let body = serde_json::to_string(&spec).expect("serialize spec");
                writeln!(w, r#"{{"name":"{name}","spec":{body}}}"#).expect("write");
                ok += 1;
            }
            Err(e) => {
                let msg = e.to_string().replace('"', "'");
                writeln!(w, r#"{{"name":"{name}","error":"{msg}"}}"#).expect("write");
                bad += 1;
            }
        }
    }
    w.flush().expect("flush");
    eprintln!("  {ok} specs built, {bad} rejected");
    ExitCode::SUCCESS
}
