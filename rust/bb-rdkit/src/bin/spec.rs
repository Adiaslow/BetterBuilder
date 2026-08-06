//! bb-spec — SMILES → MoleculeSpec JSON (on the patched RDKit, pure Rust FFI). Pipe into bb-embed.
//!
//!   bb-spec "<smiles>" > spec.json
//! (the RDKit dylib path is baked via rpath at build time, so it usually just runs.)

use std::process::ExitCode;

fn main() -> ExitCode {
    let Some(smiles) = std::env::args().nth(1) else {
        eprintln!("usage: bb-spec <smiles>");
        return ExitCode::from(2);
    };
    match bb_rdkit::build_spec(&smiles) {
        Ok(spec) => {
            println!("{}", serde_json::to_string(&spec).expect("serialize spec"));
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("bb-spec: {e}");
            ExitCode::FAILURE
        }
    }
}
