//! bb-embed — MoleculeSpec JSON (stdin or file) → conformer coords JSON (stdout).
//!
//!   bb-embed <spec.json> [mode] [seed]  > coords.json     mode = "recipe" (default) | "<n_conf>"
//!   bb-spec-native "<smiles>" | bb-embed - recipe 210185 > coords.json
//!   bb-embed <spec.json> score < coords.json  > energies.json   (validation: score EXTERNAL coords)
//!
//! "recipe" runs the two-stage core-pin recipe (conformer count from the spec); a number runs that
//! many independent embeds instead. Output: {"n_atoms": N, "conformers": [[x0,y0,z0, x1,...], ...]}
//! (n_atoms*3 per conformer, incl. H; atom order = bb-rdkit's SmilesToMol+AddHs).
//!
//! "score" reads conformer coords from stdin ({"conformers":[[...]]} or a bare [[...]]; 3 per atom)
//! and emits {"energies":[...]} — each the 3D Stage-A bounds-violation energy of that conformer on
//! the spec's (certified-identical) DistGeom bounds matrix. This is a neutral geometric ruler that
//! scores ANY ensemble (hers or native's) against the shared bounds both pipelines embed to.

use std::io::Read;

use bb_core::MoleculeSpec;
use serde::Serialize;

#[derive(Serialize)]
struct Out {
    n_atoms: usize,
    conformers: Vec<Vec<f64>>,
}

fn read_input(path: &str) -> String {
    if path == "-" {
        let mut s = String::new();
        std::io::stdin().read_to_string(&mut s).expect("read stdin");
        s
    } else {
        std::fs::read_to_string(path).expect("read file")
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).map(String::as_str).unwrap_or("-");
    let mode = args.get(2).map(String::as_str).unwrap_or("recipe");
    let seed: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(210185);

    if mode == "score" {
        // spec from the path arg (file); conformer coords JSON from stdin.
        let spec: MoleculeSpec =
            serde_json::from_str(&read_input(path)).expect("parse MoleculeSpec JSON");
        let mut cs = String::new();
        std::io::stdin().read_to_string(&mut cs).expect("read coords stdin");
        let v: serde_json::Value = serde_json::from_str(&cs).expect("parse coords JSON");
        let confs: Vec<Vec<f64>> = match v {
            serde_json::Value::Object(m) => {
                serde_json::from_value(m["conformers"].clone()).expect("coords.conformers")
            }
            other => serde_json::from_value(other).expect("bare conformer list"),
        };
        let energies: Vec<f64> = confs
            .iter()
            .map(|c| {
                bb_embed::forcefield::stage_a_energy_grad(
                    &spec,
                    c,
                    3,
                    bb_embed::forcefield::BASIN_ALL,
                )
                .0
            })
            .collect();
        println!(
            "{}",
            serde_json::to_string(&serde_json::json!({ "energies": energies })).expect("serialize")
        );
        return;
    }

    let spec: MoleculeSpec = serde_json::from_str(&read_input(path)).expect("parse MoleculeSpec JSON");

    let embedded = match mode.parse::<usize>() {
        // independent embeds
        Ok(n_conf) => bb_embed::embed(&spec, n_conf, seed).map_err(bb_embed::RecipeError::from),
        // core-pin recipe, block after block
        Err(_) => bb_embed::embed_recipe(&spec, seed)
            .map(|blocks| blocks.into_iter().flat_map(|b| b.conformers).collect()),
    };
    let confs: Vec<bb_embed::Conformer> = match embedded {
        Ok(confs) => confs,
        Err(e) => {
            eprintln!("bb-embed: {e}");
            std::process::exit(1);
        }
    };
    let out = Out {
        n_atoms: spec.n_atoms,
        conformers: confs.into_iter().map(|c| c.coords).collect(),
    };
    println!("{}", serde_json::to_string(&out).expect("serialize"));
}
