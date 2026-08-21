//! bb-aromaticity-dump — RDKit's perceived aromaticity for a whole corpus, one JSON line per molecule.
//!
//!   bb-aromaticity-dump <corpus.smi> > arom.jsonl
//!
//! Each line: {"name":..,"atoms":[0/1..],"bonds":[[begin,end,0/1]..]} — the ground truth for the
//! bb-perceive aromaticity-perception port. A SMILES RDKit rejects is emitted as {"name":..,"reject":true}.

use std::io::{BufRead, BufWriter, Write};
use std::process::ExitCode;

fn main() -> ExitCode {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: bb-aromaticity-dump <corpus.smi>");
        return ExitCode::from(2);
    };
    let Ok(f) = std::fs::File::open(&path) else {
        eprintln!("cannot open {path}");
        return ExitCode::from(2);
    };
    let stdout = std::io::stdout();
    let mut w = BufWriter::new(stdout.lock());
    for line in std::io::BufReader::new(f).lines() {
        let Ok(line) = line else { break };
        let mut it = line.split_whitespace();
        let (Some(smiles), Some(name)) = (it.next(), it.next()) else { continue };
        let v = bb_rdkit::aromatic_perception(smiles).unwrap_or_default();
        if v.first() == Some(&-1) || v.is_empty() {
            writeln!(w, r#"{{"name":"{name}","reject":true}}"#).expect("write");
            continue;
        }
        let na = v[0] as usize;
        let atoms = &v[1..1 + na];
        let bonds: Vec<[i32; 3]> = v[1 + na..].chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
        let atoms_j = serde_json::to_string(atoms).unwrap();
        let bonds_j = serde_json::to_string(&bonds).unwrap();
        writeln!(w, r#"{{"name":"{name}","atoms":{atoms_j},"bonds":{bonds_j}}}"#).expect("write");
    }
    ExitCode::SUCCESS
}
