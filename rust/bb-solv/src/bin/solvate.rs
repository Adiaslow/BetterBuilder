//! bb-solvate — MoleculeSpec JSON + conformer coords JSON → AMSOL solvation for one conformer.
//!
//!   bb-spec "<smiles>" > spec.json
//!   bb-embed spec.json recipe 210185 > confs.json
//!   bb-solvate spec.json confs.json <outdir> [conformer index, default 0]
//!
//! Writes `temp.in-wat`, `temp.in-hex`, `temp.o-wat`, `temp.o-hex` and `output.solv` into `<outdir>`.
//! The AMSOL executable comes from the configuration described in [`bb_solv::config`].

use std::path::PathBuf;
use std::process::ExitCode;

use bb_core::MoleculeSpec;
use bb_solv::config::{AmsolFile, SolvConfig};
use bb_solv::run::AmsolRunner;
use bb_solv::solv;
use bb_solv::zmatrix::to_zmatrix;
use serde::Deserialize;

#[derive(Deserialize)]
struct Conformers {
    n_atoms: usize,
    conformers: Vec<Vec<f64>>,
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: bb-solvate <spec.json> <confs.json> <outdir> [conf-index]");
        return ExitCode::from(2);
    }
    let idx: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
    let outdir = PathBuf::from(&args[3]);

    match run(&args[1], &args[2], &outdir, idx) {
        Ok(path) => {
            println!("{}", path.display());
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("bb-solvate: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(
    spec_path: &str,
    confs_path: &str,
    outdir: &PathBuf,
    idx: usize,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let spec: MoleculeSpec = serde_json::from_str(&std::fs::read_to_string(spec_path)?)?;
    let confs: Conformers = serde_json::from_str(&std::fs::read_to_string(confs_path)?)?;

    if spec.atomic_numbers.len() != spec.n_atoms {
        return Err(format!(
            "spec carries {} atomic numbers for {} atoms; rebuild it with a bb-spec that emits them",
            spec.atomic_numbers.len(),
            spec.n_atoms
        )
        .into());
    }
    if confs.n_atoms != spec.n_atoms {
        return Err(format!(
            "spec has {} atoms, conformers have {}",
            spec.n_atoms, confs.n_atoms
        )
        .into());
    }
    let coords = confs
        .conformers
        .get(idx)
        .ok_or_else(|| format!("conformer {idx} of {}", confs.conformers.len()))?;

    let zlines = to_zmatrix(&spec.atomic_numbers, coords, &spec.bonds)?;
    let cfg = SolvConfig::from_env(&AmsolFile::default())?;
    std::fs::create_dir_all(outdir)?;

    let name = outdir
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "mol".to_string());
    let run = AmsolRunner::new(cfg).solvate(outdir, &name, spec.formal_charge, &zlines)?;

    let water = solv::parse(&std::fs::read_to_string(&run.water)?)?;
    let hexadecane = solv::parse(&std::fs::read_to_string(&run.hexadecane)?)?;
    let out = outdir.join("output.solv");
    std::fs::write(&out, solv::combine(&water, &hexadecane)?.render())?;
    Ok(out)
}
