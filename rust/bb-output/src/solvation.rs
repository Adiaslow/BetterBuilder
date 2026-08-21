//! Charges from solvation — run `bb-solv` (AMSOL) on one conformation and return the
//! [`bb_solv::solv::SolvFile`] the output serializers consume. The AMSOL charges/desolvation are
//! per-atom (conformation independent), so a single solvation conformer supplies the numbers for the
//! whole ensemble, exactly as `build_ligands.py` solvates one RDKit conformer and reuses it for the
//! db2. `SolvFile` is the one `.solv` model — `bb-solv` owns the format (its `render` is byte-gated
//! against real `.solv` files), and the serializers read the struct directly, so nothing re-parses.

use bb_core::MoleculeSpec;
use bb_solv::config::{AmsolFile, SolvConfig};
use bb_solv::run::AmsolRunner;
use bb_solv::solv::SolvFile;
use std::path::Path;

/// Solvate one conformer (`coords`, a flat `x,y,z,…` list) and return its charges + desolvation.
pub fn solvate(spec: &MoleculeSpec, coords: &[f64], workdir: &Path) -> Result<SolvFile, String> {
    let zlines = bb_solv::zmatrix::to_zmatrix(&spec.atomic_numbers, coords, &spec.bonds)
        .map_err(|e| format!("z-matrix: {e}"))?;
    let cfg = SolvConfig::from_env(&AmsolFile::default()).map_err(|e| format!("amsol config: {e}"))?;
    std::fs::create_dir_all(workdir).map_err(|e| e.to_string())?;
    let name = workdir.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| "mol".into());
    let run = AmsolRunner::new(cfg)
        .solvate(workdir, &name, spec.formal_charge, &zlines)
        .map_err(|e| format!("amsol run: {e}"))?;
    let water = bb_solv::solv::parse(&std::fs::read_to_string(&run.water).map_err(|e| e.to_string())?)
        .map_err(|e| format!("parse water: {e}"))?;
    let hexadecane = bb_solv::solv::parse(&std::fs::read_to_string(&run.hexadecane).map_err(|e| e.to_string())?)
        .map_err(|e| format!("parse hexadecane: {e}"))?;
    bb_solv::solv::combine(&water, &hexadecane).map_err(|e| format!("combine: {e}"))
}
