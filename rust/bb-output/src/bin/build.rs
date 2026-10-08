//! bb-build — the RDKit-free drop-in for `build_ligands.py`: a SMILES file → the DOCK output tarball.
//!
//! ```text
//! bb-build <input.smi> <output.tar.gz> <workdir> [seed]
//! ```
//!
//! `input.smi` is one `smiles name [prot_id]` per line ([`bb_output::assemble::read_input`]). Each
//! molecule runs the full pipeline (perceive → embed → AMSOL solvate → SYBYL type + strain → mol2 +
//! db2); all members are packed into one gzip tar, each filed by its name
//! ([`bb_output::tarball::archive_dir`]). AMSOL is configured from the environment (`BB_AMSOL_EXE`,
//! `BB_AMSOL_LD_LIBRARY_PATH`). A line or molecule that fails — no name, a name that cannot name files,
//! a bad protomer id, or a pipeline error — is reported to stderr and skipped, so one bad input does
//! not sink the batch.

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: bb-build <input.smi> <output.tar.gz> <workdir> [seed]");
        return ExitCode::from(2);
    }
    let (input, output, workdir) = (&args[1], &args[2], std::path::PathBuf::from(&args[3]));
    let seed: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0x0BADC0DE);

    let text = match std::fs::read_to_string(input) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("bb-build: cannot read {input}: {e}");
            return ExitCode::from(2);
        }
    };

    let mut ligands = Vec::new();
    let (mut ok, mut fail) = (0usize, 0usize);
    for entry in bb_output::assemble::read_input(&text) {
        let bb_output::assemble::InputMolecule { smiles, name, prot_id } = match entry {
            Ok(m) => m,
            Err(e) => {
                eprintln!("bb-build: {e}");
                fail += 1;
                continue;
            }
        };
        let wd = workdir.join(format!("{name}.{prot_id}"));
        match bb_output::assemble::ligand_from_smiles(&name, prot_id, &smiles, seed, &wd) {
            Ok(lig) => {
                ligands.push(lig);
                ok += 1;
            }
            Err(e) => {
                eprintln!("bb-build: {name} failed: {e}");
                fail += 1;
            }
        }
    }

    match bb_output::write_tarball(&ligands) {
        Ok(bytes) => {
            if let Err(e) = std::fs::write(output, bytes) {
                eprintln!("bb-build: cannot write {output}: {e}");
                return ExitCode::FAILURE;
            }
        }
        Err(e) => {
            eprintln!("bb-build: tarball error: {e}");
            return ExitCode::FAILURE;
        }
    }
    eprintln!("bb-build: {ok} built, {fail} failed → {output}");
    if ok > 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
