//! Gate the mol2db2 port against the container's own `mol2db2`, per fixture in
//! `validation/db2/reference/`: `<base>.json` (input + every `hierarchy.Hierarchy` internal, dumped by
//! `dump_mol2db2.py`) and `<base>.db2` (golden output of the container's `mol2db2.py`). Each internal
//! array is checked (so a failure localizes to a stage), then the assembled db2 is required to be
//! byte-identical to golden. Runs under `cargo test -p bb-output`.

use std::path::PathBuf;

use bb_output::{build, hierarchy, mol::TypedMol};
use bb_solv::solv::{SolvAtom, SolvFile};
use serde::Deserialize;

#[derive(Deserialize)]
struct Input {
    name: String,
    #[serde(rename = "protName")]
    prot_name: String,
    smiles: String,
    longname: String,
    #[serde(rename = "atomNum")]
    atom_num: Vec<usize>,
    #[serde(rename = "atomName")]
    atom_name: Vec<String>,
    #[serde(rename = "atomType")]
    atom_type: Vec<String>,
    #[serde(rename = "dockNum")]
    dock_num: Vec<i32>,
    #[serde(rename = "colorNum")]
    color_num: Vec<i32>,
    #[serde(rename = "bondNum")]
    bond_num: Vec<usize>,
    #[serde(rename = "bondStart")]
    bond_start: Vec<usize>,
    #[serde(rename = "bondEnd")]
    bond_end: Vec<usize>,
    #[serde(rename = "bondType")]
    bond_type: Vec<String>,
    #[serde(rename = "atomXyz")]
    atom_xyz: Vec<Vec<[f64; 3]>>,
    #[serde(rename = "inputTotalStrain")]
    input_total_strain: Vec<f64>,
    #[serde(rename = "inputMaxStrain")]
    input_max_strain: Vec<f64>,
    #[serde(rename = "inputHydrogens")]
    input_hydrogens: Vec<i32>,
}

#[derive(Deserialize)]
struct SolvJson {
    #[serde(rename = "totalCharge")]
    total_charge: f64,
    #[serde(rename = "totalPolarSolv")]
    total_polar_solv: f64,
    #[serde(rename = "totalApolarSolv")]
    total_apolar_solv: f64,
    #[serde(rename = "totalSolv")]
    total_solv: f64,
    #[serde(rename = "totalSurface")]
    total_surface: f64,
    charge: Vec<f64>,
    #[serde(rename = "polarSolv")]
    polar_solv: Vec<f64>,
    #[serde(rename = "apolarSolv")]
    apolar_solv: Vec<f64>,
    solv: Vec<f64>,
    surface: Vec<f64>,
}

#[derive(Deserialize)]
struct Hier {
    #[serde(rename = "rigidStructures")]
    rigid_structures: Vec<usize>,
    #[serde(rename = "posCount")]
    pos_count: Vec<usize>,
    #[serde(rename = "numConfs")]
    num_confs: usize,
    #[serde(rename = "heavyRigidAtomNums")]
    heavy_rigid_atom_nums: Vec<usize>,
    #[serde(rename = "outAtomOrigAtom")]
    out_atom_orig_atom: Vec<usize>,
    #[serde(rename = "outAtomConfNum")]
    out_atom_conf_num: Vec<usize>,
    #[serde(rename = "brokenSets")]
    broken_sets: Vec<usize>,
}

#[derive(Deserialize)]
struct Fixture {
    input: Input,
    solv: SolvJson,
    hier: Hier,
}

fn ref_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validation/db2/reference")
}

fn load(base: &str) -> (TypedMol, SolvFile, Hier) {
    let dir = ref_dir();
    let j: Fixture = serde_json::from_str(&std::fs::read_to_string(dir.join(format!("{base}.json"))).unwrap()).unwrap();
    let mut m = TypedMol {
        name: j.input.name,
        prot_name: j.input.prot_name,
        smiles: j.input.smiles,
        longname: j.input.longname,
        atom_num: j.input.atom_num,
        atom_name: j.input.atom_name,
        atom_type: j.input.atom_type,
        atom_bonds: Vec::new(),
        bond_num: j.input.bond_num,
        bond_start: j.input.bond_start,
        bond_end: j.input.bond_end,
        bond_type: j.input.bond_type,
        atom_xyz: j.input.atom_xyz,
        input_total_strain: j.input.input_total_strain,
        input_max_strain: j.input.input_max_strain,
        input_hydrogens: j.input.input_hydrogens,
    };
    m.build_adjacency();
    // Map the fixture's `.solv` fields into the one shared model (bb-solv::SolvFile): the per-atom
    // columns are charge / diff_pol(polarSolv) / area(surface) / diff_apol(apolarSolv) /
    // diff_atomic_solv(solv), matching the `.solv` layout bb-solv renders.
    let atoms = (0..j.solv.charge.len())
        .map(|i| SolvAtom {
            charge: j.solv.charge[i],
            diff_pol: j.solv.polar_solv[i],
            area: j.solv.surface[i],
            diff_apol: j.solv.apolar_solv[i],
            diff_atomic_solv: j.solv.solv[i],
        })
        .collect();
    let s = SolvFile {
        name: String::new(),
        charge: j.solv.total_charge,
        area: j.solv.total_surface,
        tot_diff_pol: j.solv.total_polar_solv,
        tot_diff_apol: j.solv.total_apolar_solv,
        tot_diff_pol_plus_apol: j.solv.total_solv,
        atoms,
    };
    // typing gate: our dockNum/colorNum must equal the container's
    for (i, &num) in m.atom_num.iter().enumerate() {
        assert_eq!(bb_output::sybyl2dock::convert_atom(&m, num).unwrap(), j.input.dock_num[i], "{base} dockNum atom {i}");
        assert_eq!(bb_output::color::convert_atom(&m, num), j.input.color_num[i], "{base} colorNum atom {i}");
    }
    (m, s, j.hier)
}

const FIXTURES: &[&str] = &[
    // single-conformation (the pipeline's solvation output.mol2)
    "benzoic.0", "mc0001.0", "mc0002.0", "mc0003.0", "mc0004.0", "mc0005.0", "mc0006.0", "mc0007.0", "mc0008.0",
    // multi-conformation (synthesized to exercise the buckets2 spatial-clustering path)
    "benzoic.0_mc", "mc0001.0_mc", "mc0003.0_mc",
    // non-contiguous movers (fragments the rigid core; adds clash-broken multi-conf sets)
    "mc0001.0_alt",
];

#[test]
fn hierarchy_internals_match() {
    for base in FIXTURES {
        let (m, _s, h) = load(base);
        let built = hierarchy::build(&m, bb_output::DISTTOL);
        assert_eq!(built.rigid_structures, h.rigid_structures, "{base} rigidStructures");
        assert_eq!(built.pos_count, h.pos_count, "{base} posCount");
        assert_eq!(built.num_confs, h.num_confs, "{base} numConfs");
        assert_eq!(built.heavy_rigid_atom_nums, h.heavy_rigid_atom_nums, "{base} heavyRigidAtomNums");
        assert_eq!(built.out_atom_orig_atom, h.out_atom_orig_atom, "{base} outAtomOrigAtom");
        assert_eq!(built.out_atom_conf_num, h.out_atom_conf_num, "{base} outAtomConfNum");
        assert_eq!(built.broken_sets, h.broken_sets, "{base} brokenSets");
    }
}

#[test]
fn tarball_members() {
    use bb_output::tarball::{full_name, ChargeLetter};
    use bb_output::MoleculeName;
    // charge char: N neutral, M -1, O +1 (chr(78 + charge))
    let name = |s: &str| MoleculeName::new(s).unwrap();
    let letter = |c: i32| ChargeLetter::new(c).unwrap();
    assert_eq!(full_name(&name("benzoic"), 0, letter(0)), "benzoic.0.N");
    assert_eq!(full_name(&name("x"), 0, letter(-1)), "x.0.M");
    assert_eq!(full_name(&name("x"), 2, letter(1)), "x.2.O");

    // assemble a real ligand's two artifacts and verify member names + byte-exact contents survive
    let (m, s, _) = load("benzoic.0");
    let charges: Vec<f64> = s.atoms.iter().map(|a| a.charge).collect();
    let mol2 = bb_output::write_mol2(&m, &charges, 0);
    let db2 = bb_db2::write_entry(&build(&m, &s).unwrap());
    let bytes = bb_output::write_tarball(&[bb_output::Ligand {
        name: name("benzoic"),
        prot_id: 0,
        charge: letter(s.charge.round() as i32),
        mol2: mol2.clone(),
        db2: db2.clone(),
    }])
    .unwrap();

    // extract and check — the real pipeline names these ./benzoic.0.N.{mol2,db2.gz}
    let dec = flate2::read::GzDecoder::new(&bytes[..]);
    let mut ar = tar::Archive::new(dec);
    let mut found = std::collections::HashMap::new();
    for e in ar.entries().unwrap() {
        let mut e = e.unwrap();
        let name = e.path().unwrap().to_string_lossy().into_owned();
        let mut content = String::new();
        std::io::Read::read_to_string(&mut e, &mut content).unwrap();
        found.insert(name, content);
    }
    // the driver writes these as `./benzoic.0.N.{ext}`; the leading `./` normalizes away on
    // extraction (POSIX), so the extracted files are identical — we check the normalized names.
    assert_eq!(found.get("benzoic.0.N.mol2"), Some(&mol2), "mol2 member");
    assert_eq!(found.get("benzoic.0.N.db2.gz"), Some(&db2), "db2 member");
    assert_eq!(found.len(), 2);
}

#[test]
fn mol2_byte_identical() {
    // the tarball's mol2 member — our writer, fed the reference atoms + solv charges, must reproduce
    // the container's `output.mol2` (written by `mol2amsol.write_mol2`) byte-for-byte.
    let dir = ref_dir();
    let mut fails = Vec::new();
    for base in FIXTURES.iter().filter(|b| !b.contains('_')) {
        let (m, s, _) = load(base);
        let charges: Vec<f64> = s.atoms.iter().map(|a| a.charge).collect();
        let ours = bb_output::write_mol2(&m, &charges, 0);
        let golden = std::fs::read_to_string(dir.join(format!("{base}.omol2"))).unwrap();
        if ours != golden {
            let i = ours.lines().zip(golden.lines()).position(|(a, b)| a != b);
            fails.push(match i {
                Some(i) => format!("{base}: line {}:\n    ours  : {:?}\n    golden: {:?}", i + 1, ours.lines().nth(i).unwrap_or("<EOF>"), golden.lines().nth(i).unwrap_or("<EOF>")),
                None => format!("{base}: length {} vs {}", ours.lines().count(), golden.lines().count()),
            });
        }
    }
    assert!(fails.is_empty(), "mol2 mismatches:\n{}", fails.join("\n"));
}

#[test]
fn db2_byte_identical() {
    let dir = ref_dir();
    let mut fails = Vec::new();
    for base in FIXTURES {
        let (m, s, _) = load(base);
        let entry = build(&m, &s).unwrap();
        let ours = bb_db2::write_entry(&entry);
        let golden = std::fs::read_to_string(dir.join(format!("{base}.db2"))).unwrap();
        if ours != golden {
            let n = ours.lines().zip(golden.lines()).position(|(a, b)| a != b);
            match n {
                Some(i) => fails.push(format!(
                    "{base}: first diff line {}:\n    ours  : {:?}\n    golden: {:?}",
                    i + 1,
                    ours.lines().nth(i).unwrap_or("<EOF>"),
                    golden.lines().nth(i).unwrap_or("<EOF>")
                )),
                None => fails.push(format!("{base}: differ in length ({} vs {} lines)", ours.lines().count(), golden.lines().count())),
            }
        }
    }
    assert!(fails.is_empty(), "db2 mismatches:\n{}", fails.join("\n"));
}

