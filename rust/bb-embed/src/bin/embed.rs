//! bb-embed — MoleculeSpec JSON (stdin or file) → candidate conformer coords JSON (stdout).
//!
//!   bb-embed <spec.json> [mode] [seed]  > coords.json     mode = "recipe" (default) | "<n_conf>"
//!   bb-spec "<smiles>" | bb-embed - recipe 210185 > coords.json
//!
//! "recipe" runs the faithful two-stage core-pin recipe (count from the spec); a number runs that
//! many independent embeds instead. Output: {"n_atoms": N, "conformers": [[x0,y0,z0, x1,...], ...]}
//! (n_atoms*3 per conformer, incl. H; atom order = bb-rdkit's SmilesToMol+AddHs).

use std::io::Read;

use bb_core::MoleculeSpec;
use serde::Serialize;

#[derive(Serialize)]
struct Out {
    n_atoms: usize,
    conformers: Vec<Vec<f64>>,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).map(String::as_str).unwrap_or("-");
    let mode = args.get(2).map(String::as_str).unwrap_or("recipe");
    let seed: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(210185);

    let text = if path == "-" {
        let mut s = String::new();
        std::io::stdin().read_to_string(&mut s).expect("read stdin");
        s
    } else {
        std::fs::read_to_string(path).expect("read spec file")
    };
    let spec: MoleculeSpec = serde_json::from_str(&text).expect("parse MoleculeSpec JSON");

    let confs = match mode.parse::<usize>() {
        Ok(n_conf) => bb_embed::embed(&spec, n_conf, seed), // independent embeds
        Err(_) => bb_embed::embed_recipe(&spec, seed),      // faithful core-pin recipe
    };
    let out = Out {
        n_atoms: spec.n_atoms,
        conformers: confs.into_iter().map(|c| c.coords).collect(),
    };
    println!("{}", serde_json::to_string(&out).expect("serialize"));
}
