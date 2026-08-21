//! Gate the native SYBYL typer against `rdkit_mol2_writer`'s own output (`validation/db2/reference/3d/`,
//! the `3d/name.mol2` files, in RDKit atom order). For each molecule we perceive its SMILES with
//! bb-perceive — whose atom order is gated identical to RDKit's — assert the element sequence matches
//! the reference (so the per-atom comparison is valid), then require every SYBYL type to match.

use std::path::PathBuf;

fn ref_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validation/db2/reference/3d")
}

fn elem(z: u8) -> &'static str {
    match z {
        1 => "H", 6 => "C", 7 => "N", 8 => "O", 9 => "F", 15 => "P", 16 => "S", 17 => "Cl", 35 => "Br", 53 => "I", _ => "*",
    }
}

/// Parse a `3d/name.mol2`'s ATOM section: returns `(element, sybyl_type)` per atom, in file order.
fn parse_3d(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut in_atoms = false;
    for line in text.lines() {
        if line.starts_with("@<TRIPOS>ATOM") {
            in_atoms = true;
            continue;
        }
        if line.starts_with("@<TRIPOS>") {
            in_atoms = false;
            continue;
        }
        if in_atoms && !line.trim().is_empty() {
            let t: Vec<&str> = line.split_whitespace().collect();
            let sybyl = t[5].to_string();
            let element = sybyl.split('.').next().unwrap().to_string();
            out.push((element, sybyl));
        }
    }
    out
}

/// Parse the ATOM names and BOND records `(a, b, type)` (1-based, sorted pair) from a 3d mol2.
fn parse_3d_names_bonds(text: &str) -> (Vec<String>, std::collections::HashSet<(usize, usize, String)>) {
    let mut names = Vec::new();
    let mut bonds = std::collections::HashSet::new();
    let mut sect = "";
    for line in text.lines() {
        if line.starts_with("@<TRIPOS>ATOM") {
            sect = "atom";
            continue;
        } else if line.starts_with("@<TRIPOS>BOND") {
            sect = "bond";
            continue;
        } else if line.starts_with("@<TRIPOS>") {
            sect = "";
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        if sect == "atom" && t.len() >= 2 {
            names.push(t[1].to_string());
        } else if sect == "bond" && t.len() >= 4 {
            let (a, b): (usize, usize) = (t[1].parse().unwrap(), t[2].parse().unwrap());
            bonds.insert((a.min(b), a.max(b), t[3].to_string()));
        }
    }
    (names, bonds)
}

#[test]
fn bond_types_and_names_match_rdkit() {
    let dir = ref_dir();
    let mut fails = Vec::new();
    for line in std::fs::read_to_string(dir.join("smiles.txt")).unwrap().lines() {
        let mut it = line.split_whitespace();
        let (smiles, name) = (it.next().unwrap(), it.next().unwrap());
        let (ref_names, ref_bonds) = parse_3d_names_bonds(&std::fs::read_to_string(dir.join(format!("{name}.mol2"))).unwrap());
        let p = bb_perceive::smarts_match::perceive(smiles).expect("perceive");
        // atom names (per-atom, order matches RDKit)
        let names = bb_output::sybyl::atom_names(&p);
        if names != ref_names {
            fails.push(format!("{name}: atom names differ"));
            continue;
        }
        // bond types as a set of (a, b, type) — bond ordering may differ, atom indices match
        let (types, amide) = bb_output::sybyl::assign_sybyl_full(&p);
        let bt = bb_output::sybyl::bond_types(&p, &types, &amide);
        let ours: std::collections::HashSet<(usize, usize, String)> = p
            .bonds
            .iter()
            .zip(&bt)
            .map(|(&(a, b), ty)| (a.min(b) + 1, a.max(b) + 1, ty.clone()))
            .collect();
        if ours != ref_bonds {
            let extra: Vec<_> = ours.difference(&ref_bonds).take(4).collect();
            let missing: Vec<_> = ref_bonds.difference(&ours).take(4).collect();
            fails.push(format!("{name}: bonds differ. ours-only {extra:?}; rdkit-only {missing:?}"));
        }
    }
    assert!(fails.is_empty(), "bond/name mismatches:\n{}", fails.join("\n"));
}

#[test]
fn sybyl_types_match_rdkit() {
    let dir = ref_dir();
    let smiles_txt = std::fs::read_to_string(dir.join("smiles.txt")).unwrap();
    let mut total_atoms = 0;
    let mut fails = Vec::new();
    for line in smiles_txt.lines() {
        let mut it = line.split_whitespace();
        let (smiles, name) = (it.next().unwrap(), it.next().unwrap());
        let ref3d = parse_3d(&std::fs::read_to_string(dir.join(format!("{name}.mol2"))).unwrap());
        let p = bb_perceive::smarts_match::perceive(smiles).expect("perceive");

        // atom-order sanity: element sequence must match (else per-atom comparison is meaningless)
        if p.atomic_numbers.len() != ref3d.len() {
            fails.push(format!("{name}: atom count {} != {}", p.atomic_numbers.len(), ref3d.len()));
            continue;
        }
        let order_ok = p.atomic_numbers.iter().zip(&ref3d).all(|(&z, (e, _))| elem(z) == e);
        if !order_ok {
            fails.push(format!("{name}: element order differs from RDKit 3d mol2"));
            continue;
        }

        let ours = bb_output::sybyl::assign_sybyl(&p);
        for (i, (o, (_, want))) in ours.iter().zip(&ref3d).enumerate() {
            total_atoms += 1;
            // exact SYBYL-string match — native Kekulization (bb_perceive::kekulize) reproduces
            // RDKit's aro6 without any heuristic, so every atom's type matches rdkit_mol2_writer.
            if o != want {
                fails.push(format!("{name} atom {i} ({}): ours {o:?} != rdkit {want:?}", ref3d[i].0));
            }
        }
    }
    println!("{total_atoms} atoms; all SYBYL types exact vs rdkit_mol2_writer");
    assert!(fails.is_empty(), "SYBYL mismatches:\n{}", fails.join("\n"));
}
