//! bb-spec — SMILES → MoleculeSpec JSON on stdout. Exit 2 on missing argument, 1 on a SMILES the
//! patched RDKit cannot parse.
//!
//!   bb-spec "<smiles>" > spec.json

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
