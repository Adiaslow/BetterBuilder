//! bb-spec-native — SMILES → MoleculeSpec JSON, built entirely from the pure-Rust ports with no
//! RDKit. Mirrors `bb-spec`'s CLI so it drops into the same pipelines. Exit 2 on missing argument,
//! 1 on a SMILES that cannot be built.

use std::process::ExitCode;

fn main() -> ExitCode {
    let Some(smiles) = std::env::args().nth(1) else {
        eprintln!("usage: bb-spec-native <smiles>");
        return ExitCode::from(2);
    };
    match bb_spec::build_native(&smiles) {
        Ok(spec) => {
            println!("{}", serde_json::to_string(&spec).expect("serialize spec"));
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("bb-spec-native: {e}");
            ExitCode::FAILURE
        }
    }
}
