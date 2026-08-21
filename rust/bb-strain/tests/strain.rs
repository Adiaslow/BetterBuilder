//! Gate the torsion-strain port against the real `Torsion_Strain` (`validation/db2/reference/`):
//! `strain_ref.json` holds the container's total/max strain for each fixture's `3d/*.mol2`. We perceive
//! the SMILES (our atom order == RDKit's == the 3d mol2's), read that mol2's coordinates, compute
//! strain, and require it to match the oracle. Same geometry in → same strain out.

use std::path::PathBuf;

fn ref_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validation/db2/reference")
}

/// Coordinates from a 3d mol2's ATOM section, in file (= perception) order.
fn parse_coords(text: &str) -> Vec<[f64; 3]> {
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
            out.push([t[2].parse().unwrap(), t[3].parse().unwrap(), t[4].parse().unwrap()]);
        }
    }
    out
}

#[test]
fn strain_matches_oracle() {
    let dir = ref_dir();
    let refs: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("strain_ref.json")).unwrap()).unwrap();
    let smiles_txt = std::fs::read_to_string(dir.join("3d/smiles.txt")).unwrap();
    let mut fails = Vec::new();
    let mut n = 0;
    for line in smiles_txt.lines() {
        let mut it = line.split_whitespace();
        let (smiles, name) = (it.next().unwrap(), it.next().unwrap());
        let want = &refs[name];
        if !want["total"].is_number() {
            continue; // oracle errored on this one
        }
        let (want_total, want_max) = (want["total"].as_f64().unwrap(), want["max"].as_f64().unwrap());

        let p = bb_perceive::smarts_match::perceive(smiles).expect("perceive");
        let coords = parse_coords(&std::fs::read_to_string(dir.join(format!("3d/{name}.mol2"))).unwrap());
        assert_eq!(coords.len(), p.atomic_numbers.len(), "{name}: coord/atom count");
        let hits = bb_strain::match_library(&p);
        let (total, max) = bb_strain::strain(&hits, &coords);

        n += 1;
        // strain energies are small (kcal/mol); require tight agreement
        if (total - want_total).abs() > 1e-3 || (max - want_max).abs() > 1e-3 {
            fails.push(format!("{name}: ours (total {total:.4}, max {max:.4}) != oracle (total {want_total:.4}, max {want_max:.4})"));
        }
    }
    assert!(n > 0, "no fixtures checked");
    assert!(fails.is_empty(), "{n} checked; strain mismatches:\n{}", fails.join("\n"));
}
