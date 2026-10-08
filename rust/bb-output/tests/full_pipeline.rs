//! Full pipeline, RDKit-free, with real AMSOL charges: SMILES → perceive → conformers → solvate
//! (AMSOL) → SYBYL type + strain → `TypedMol` → mol2 + db2 → tarball. Needs the AMSOL toolchain, so it
//! is ignored by default — reported as skipped, never as passed — and runs with `--ignored` where AMSOL
//! exists (`rust/gates.sh dbgate`); run that way without AMSOL, it fails and says so.

use std::path::Path;

/// Point `bb-solv` at AMSOL: `BB_AMSOL_EXE` if set, else the toolchain AMSOL on Wynton. `false` when
/// neither exists.
fn setup_amsol() -> bool {
    if std::env::var_os("BB_AMSOL_EXE").is_some() {
        return true;
    }
    let exe = Path::new("/nfs/home/amurray2/toolchain/amsol/amsol7.1");
    let lib = "/nfs/home/amurray2/toolchain/amsol/lib";
    if !exe.exists() {
        return false;
    }
    std::env::set_var("BB_AMSOL_EXE", exe);
    std::env::set_var("BB_AMSOL_LD_LIBRARY_PATH", lib);
    true
}

#[test]
#[ignore = "needs AMSOL (BB_AMSOL_EXE, or the Wynton toolchain); run with --ignored, as rust/gates.sh dbgate does"]
fn smiles_to_tarball_with_amsol_charges() {
    assert!(setup_amsol(), "AMSOL not found: set BB_AMSOL_EXE (and BB_AMSOL_LD_LIBRARY_PATH)");
    // `build_tarball_from_smiles` names the member `{name}.{prot_id}.{C}` with prot_id=0, so `name`
    // must be the BARE molecule name (not "benzoic.0" — that would double the prot id to benzoic.0.0.N).
    let workdir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("benzoic");
    let bytes = bb_output::assemble::build_tarball_from_smiles(&bb_output::MoleculeName::new("benzoic").unwrap(), "OC(=O)c1ccccc1", 0xBEEF, &workdir)
        .expect("build tarball with AMSOL charges");

    // extract and check both members exist and the mol2 carries real (non-zero) charges
    let dec = flate2::read::GzDecoder::new(&bytes[..]);
    let mut ar = tar::Archive::new(dec);
    let mut members = std::collections::HashMap::new();
    for e in ar.entries().unwrap() {
        let mut e = e.unwrap();
        let name = e.path().unwrap().to_string_lossy().into_owned();
        let mut c = String::new();
        std::io::Read::read_to_string(&mut e, &mut c).unwrap();
        members.insert(name, c);
    }
    let mol2 = members.get("benzoic.0.N.mol2").expect("mol2 member");
    assert!(members.contains_key("benzoic.0.N.db2.gz"), "db2 member");

    // the ATOM charge column (last field) must not be all zeros — AMSOL produced real charges
    let any_nonzero = mol2
        .lines()
        .skip_while(|l| !l.starts_with("@<TRIPOS>ATOM"))
        .take_while(|l| !l.starts_with("@<TRIPOS>BOND"))
        .filter_map(|l| l.split_whitespace().last())
        .filter_map(|q| q.parse::<f64>().ok())
        .any(|q| q.abs() > 1e-6);
    assert!(any_nonzero, "AMSOL charges should be non-zero in the mol2");

    // and the db2 round-trips — ONE entry per core-seed block, as build_ligands.py concatenates its
    // per-block `db2ins` (benzoic: core_seeds=20 → 20 entries). A single entry over all confs would be
    // the pre-per-block bug.
    let db2 = members.get("benzoic.0.N.db2.gz").unwrap();
    assert_eq!(bb_db2::parse(db2).unwrap().len(), 20, "one db2 entry per core-seed block");
}
