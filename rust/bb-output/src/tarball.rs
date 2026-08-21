//! Assemble the pipeline's output tarball — a faithful reproduction of `build_ligands.py`'s
//! `output.tar.gz`. Each built ligand contributes two members, named exactly as the driver's
//! `tar_name(get_zinc_directory_hash(name), fullname, ext)`: a `.mol2` (the charge-bearing mol2) and a
//! `.db2.gz` (the db2 text — the driver stores it *uncompressed* under a `.gz` name, so we match that).
//! `fullname` is `name.prot_id.C` where `C = chr(78 + formal_charge)` (N neutral, M −1, O +1, …). The
//! archive is a gzipped tar. The driver's own tar is left truncated by a close-ordering bug; ours is
//! finalized properly — the one intentional deviation, a fix not a divergence.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;

/// One ligand's output: identity + the two already-built artifact texts.
pub struct Ligand {
    pub name: String,
    pub prot_id: i32,
    /// net formal charge (from the `.solv`/db2 total charge, rounded)
    pub formal_charge: i32,
    pub mol2: String,
    pub db2: String,
}

/// `mol_fullname` = `name.prot_id.C`, `C = chr(78 + formal_charge)`.
pub fn full_name(name: &str, prot_id: i32, formal_charge: i32) -> String {
    let c = char::from((78 + formal_charge) as u8);
    format!("{name}.{prot_id}.{c}")
}

/// The tarball directory prefix for a molecule name. ZINC ids hash into a tranche path; everything
/// else (our macrocycles) goes to `.` — reproducing `get_zinc_directory_hash`. (The ZINC tranche
/// branch is not needed for the current corpus and is left unimplemented rather than guessed.)
pub fn zinc_hash(name: &str) -> &'static str {
    if name.starts_with("ZINC") {
        unimplemented!("ZINC tranche hashing not needed for the macrocycle corpus");
    }
    "."
}

/// `tar_name(hash, fullname, ext)` = `hash + '/' + fullname + '.' + ext` (or just `fullname.ext` when
/// hash is empty).
fn tar_name(hash: &str, fullname: &str, ext: &str) -> String {
    if hash.is_empty() {
        format!("{fullname}.{ext}")
    } else {
        format!("{hash}/{fullname}.{ext}")
    }
}

fn add(builder: &mut tar::Builder<impl Write>, name: &str, data: &[u8]) -> std::io::Result<()> {
    let mut header = tar::Header::new_gnu();
    header.set_size(data.len() as u64);
    header.set_mode(0o644);
    header.set_mtime(0); // deterministic: the tar carries no wall-clock time
    header.set_cksum();
    builder.append_data(&mut header, name, data)
}

/// Build the gzipped output tarball for a set of ligands, returning its bytes. Members are emitted in
/// the given order; for each ligand a `.mol2` then a `.db2.gz` (matching the driver's per-ligand order).
pub fn write_tarball(ligands: &[Ligand]) -> std::io::Result<Vec<u8>> {
    let gz = GzEncoder::new(Vec::new(), Compression::default());
    let mut builder = tar::Builder::new(gz);
    for lig in ligands {
        let hash = zinc_hash(&lig.name);
        let fullname = full_name(&lig.name, lig.prot_id, lig.formal_charge);
        add(&mut builder, &tar_name(hash, &fullname, "mol2"), lig.mol2.as_bytes())?;
        add(&mut builder, &tar_name(hash, &fullname, "db2.gz"), lig.db2.as_bytes())?;
    }
    let gz = builder.into_inner()?;
    gz.finish()
}
