//! End-to-end: SMILES → BetterBuilder's own pipeline → the DOCK output tarball, fully RDKit-free.
//! Perceive → SYBYL type → embed our own conformers (ETKDG core-pin) → assemble a `TypedMol` →
//! serialize both members (mol2 + db2) → tar. This is a structural check (our geometry differs from
//! Divya's, as established for solvation): the tarball is well-formed, its two members are self-
//! consistent, and the db2 round-trips through `bb-db2`. The *components* are byte-gated elsewhere
//! (SYBYL vs rdkit_mol2_writer, db2/mol2 format vs the container's writers).

use bb_solv::solv::{SolvAtom, SolvFile};

/// A zero-charge solv of the given atom count (production charges come from `bb-solv`/AMSOL; geometry,
/// typing and structure are what this test exercises).
fn zero_solv(n: usize) -> SolvFile {
    SolvFile {
        name: "fake".into(),
        charge: 0.0,
        area: 0.0,
        tot_diff_pol: 0.0,
        tot_diff_apol: 0.0,
        tot_diff_pol_plus_apol: 0.0,
        atoms: vec![SolvAtom { charge: 0.0, diff_pol: 0.0, area: 0.0, diff_apol: 0.0, diff_atomic_solv: 0.0 }; n],
    }
}

#[test]
fn smiles_to_tarball() {
    // a real macrocycle from the corpus + a small molecule
    for (name, smiles) in [
        ("benzoic.0", "OC(=O)c1ccccc1"),
        ("mc0005.0", "COC(CS(=O)(=O)N[C@@H]1CC=CCC(C)C(=O)N2C[C@@H](NC(=O)C1)[C@H](CC1CC1)C2)C1CCC1"),
    ] {
        let mol = bb_output::assemble::from_smiles(name, smiles, 0xC0FFEE).expect("from_smiles");
        let n = mol.atom_num.len();
        assert!(!mol.atom_xyz.is_empty(), "{name}: got conformers");
        assert_eq!(mol.atom_type.len(), n);
        let solv = zero_solv(n);

        // both members from the one TypedMol
        let mol2 = bb_output::write_mol2(&mol, &vec![0.0; n], 0);
        let db2entry = bb_output::build(&mol, &solv).expect("build db2");
        let db2 = bb_db2::write_entry(&db2entry);

        // db2 must round-trip (well-formed, self-consistent)
        let parsed = bb_db2::parse(&db2).expect("db2 parses");
        assert_eq!(parsed.len(), 1, "{name}: one db2 entry");
        assert_eq!(parsed[0].atoms.len(), n, "{name}: db2 atom count");
        assert_eq!(bb_db2::write_entry(&parsed[0]), db2, "{name}: db2 round-trips byte-identically");

        // mol2 is well-formed
        assert!(mol2.contains("@<TRIPOS>ATOM") && mol2.contains("@<TRIPOS>BOND"), "{name}: mol2 sections");
        assert_eq!(mol2.lines().filter(|l| l.starts_with("@<TRIPOS>ATOM")).count(), 1);

        // Internal consistency — the invariant Divya's pipeline VIOLATES (her mol2 member is
        // RDKit-typed, her db2 member OpenEye-typed, so they disagree on ~1 atom each). Both of our
        // members serialize from one TypedMol, so the db2's dock vdw type must be exactly what the
        // mol2 member's SYBYL string maps to via sybyl2dock — per atom, in the same order.
        let mol2_sybyl: Vec<&str> = mol2
            .lines()
            .skip_while(|l| !l.starts_with("@<TRIPOS>ATOM"))
            .skip(1)
            .take_while(|l| !l.starts_with("@<TRIPOS>"))
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.split_whitespace().nth(5).expect("mol2 SYBYL column"))
            .collect();
        assert_eq!(mol2_sybyl.len(), n, "{name}: mol2 lists every atom");
        // range loop is intentional: each iteration touches three parallel arrays plus the 1-based
        // atom number, so indexing reads clearer than zipping.
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            assert_eq!(mol2_sybyl[i], mol.atom_type[i], "{name}: mol2 SYBYL atom {i}");
            let expected = bb_output::sybyl2dock::convert_atom(&mol, i + 1).expect("dock type");
            assert_eq!(parsed[0].atoms[i].vdwtype, expected,
                "{name}: db2 dock type of atom {i} ({}) disagrees with its mol2 SYBYL type", mol.atom_type[i]);
        }

        // assemble the tarball and verify both members extract
        let bytes = bb_output::write_tarball(&[bb_output::Ligand {
            name: bb_output::MoleculeName::new(name.split('.').next().unwrap()).unwrap(),
            prot_id: 0,
            charge: bb_output::tarball::ChargeLetter::new(0).unwrap(),
            mol2: mol2.clone(),
            db2: db2.clone(),
        }])
        .expect("tarball");
        let dec = flate2::read::GzDecoder::new(&bytes[..]);
        let mut ar = tar::Archive::new(dec);
        let members: Vec<String> = ar.entries().unwrap().map(|e| e.unwrap().path().unwrap().to_string_lossy().into_owned()).collect();
        let base = name.split('.').next().unwrap();
        assert!(members.contains(&format!("{base}.0.N.mol2")), "{name}: mol2 member ({members:?})");
        assert!(members.contains(&format!("{base}.0.N.db2.gz")), "{name}: db2 member");
    }
}
