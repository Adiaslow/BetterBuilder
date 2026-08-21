//! bb-emit-inputs — writer-gate helper: SMILES → the (mol2, solv) that native's writer consumed, plus
//! native's own db2, so the corpus-scale writer gate can feed the SAME (mol2, solv) to the container's
//! `mol2db2.py` and byte-diff its db2 against native's. Stage-isolated: identical input to both writers.
//!
//!   bb-emit-inputs <smiles> <name> <outdir> [seed]
//!   → <outdir>/<name>.mol2   (native charge-bearing mol2, `-m` for mol2db2.py)
//!     <outdir>/<name>.solv   (native rendered .solv, `-s` for mol2db2.py)
//!     <outdir>/<name>.native.db2  (native writer output, the comparison target)
//! Exit 0 on success, 1 on any pipeline error (skipped molecules are the caller's to note).
use std::path::{Path, PathBuf};
use std::process::ExitCode;

fn run(smiles: &str, name: &str, outdir: &Path, seed: u64) -> Result<(), String> {
    let (spec, p) = bb_spec::build_native_perceived(smiles)?;
    let raw = bb_embed::embed_recipe(&spec, seed);
    if raw.is_empty() {
        return Err("no conformers embedded".into());
    }
    let confs = bb_output::assemble::confs_from_embed(&raw, spec.n_atoms);
    let workdir = outdir.join(format!("{name}.amsol"));
    let mut solv = bb_output::solvation::solvate(&spec, &raw[0].coords, &workdir)?;
    if solv.atoms.len() != spec.n_atoms {
        return Err(format!("solv atoms {} != {}", solv.atoms.len(), spec.n_atoms));
    }
    // Round solv + coords to the FILE precisions her mol2db2 parses (rendered .solv / mol2), so native's
    // db2 is built from the SAME values — the single source is `assemble`, shared with production's
    // `build_blocks_db2` so the gate can never validate a differently-rounded artifact than ships.
    bb_output::assemble::round_solv_to_file_precision(&mut solv);
    let full = format!("{name}.0");
    let mut mol = bb_output::assemble::assemble(&full, smiles, &p, confs);
    bb_output::assemble::round_coords_to_file_precision(&mut mol);
    let charges: Vec<f64> = solv.atoms.iter().map(|a| a.charge).collect();
    std::fs::write(outdir.join(format!("{name}.solv")), solv.render()).map_err(|e| e.to_string())?;

    // Per block (same chunking as production), emit each block's mol2 (her mol2db2 input) + native db2.
    bb_output::assemble::for_each_block(&mol, spec.sidechain_confs as usize, |b, bmol| {
        let bdb2 = bb_db2::write_entry(&bb_output::build(bmol, &solv)?);
        let bmol2 = bb_output::write_mol2_all(bmol, &charges);
        std::fs::write(outdir.join(format!("{name}.b{b}.mol2")), bmol2).map_err(|e| e.to_string())?;
        std::fs::write(outdir.join(format!("{name}.b{b}.native.db2")), bdb2).map_err(|e| e.to_string())?;
        Ok::<(), String>(())
    })?;
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: bb-emit-inputs <smiles> <name> <outdir> [seed]");
        return ExitCode::from(2);
    }
    let outdir = PathBuf::from(&args[3]);
    let seed: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(210185);
    if let Err(e) = std::fs::create_dir_all(&outdir) {
        eprintln!("bb-emit-inputs: mkdir {}: {e}", outdir.display());
        return ExitCode::from(1);
    }
    match run(&args[1], &args[2], &outdir, seed) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("bb-emit-inputs: {} failed: {e}", args[2]);
            ExitCode::from(1)
        }
    }
}
