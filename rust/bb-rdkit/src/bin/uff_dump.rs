//! bb-uff-dump — RDKit UFF oracle capture.
//!
//!   bb-uff-dump table               → RDKit's UFF parameter table verbatim on stdout
//!   bb-uff-dump labels <corpus.smi> → per-atom UFF atom-type labels as JSONL on stdout
//!
//! `table` is the exact `defaultParamData` string the UFF typer parses, for vendoring. `labels`
//! emits one object per corpus line — {"name": <id>, "labels": ["C_3", "N_R", ...]} in post-AddHs
//! order — the oracle the Rust getAtomLabel port is gated against. A SMILES RDKit cannot parse is
//! {"name": <id>, "error": "unparseable"} rather than dropped, so the count matches the corpus.

use std::io::{BufRead, BufWriter, Write};
use std::process::ExitCode;

fn main() -> ExitCode {
    let mode = std::env::args().nth(1).unwrap_or_default();
    match mode.as_str() {
        "table" => {
            let t = bb_rdkit::uff_param_table().expect("dump UFF param table");
            print!("{t}");
            ExitCode::SUCCESS
        }
        "chiraltags" => {
            let Some(path) = std::env::args().nth(2) else {
                eprintln!("usage: bb-uff-dump chiraltags <corpus.smi>");
                return ExitCode::from(2);
            };
            let f = std::fs::File::open(&path).expect("open corpus");
            let stdout = std::io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            for line in std::io::BufReader::new(f).lines() {
                let Ok(line) = line else { break };
                let mut it = line.split_whitespace();
                let (Some(smiles), Some(name)) = (it.next(), it.next()) else { continue };
                match bb_rdkit::chiral_tags(smiles) {
                    Ok(v) if !v.is_empty() => {
                        let q: Vec<String> = v.iter().map(|x| x.to_string()).collect();
                        writeln!(w, r#"{{"name":"{name}","tags":[{}]}}"#, q.join(",")).ok();
                    }
                    _ => { writeln!(w, r#"{{"name":"{name}","error":"unparseable"}}"#).ok(); }
                }
            }
            ExitCode::SUCCESS
        }
        "bondstereo" => {
            let Some(path) = std::env::args().nth(2) else {
                eprintln!("usage: bb-uff-dump bondstereo <corpus.smi>");
                return ExitCode::from(2);
            };
            let f = std::fs::File::open(&path).expect("open corpus");
            let stdout = std::io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            for line in std::io::BufReader::new(f).lines() {
                let Ok(line) = line else { break };
                let mut it = line.split_whitespace();
                let (Some(smiles), Some(name)) = (it.next(), it.next()) else { continue };
                match bb_rdkit::bond_stereo(smiles) {
                    Ok(v) if !v.is_empty() => {
                        let q: Vec<String> = v.iter().map(|x| x.to_string()).collect();
                        writeln!(w, r#"{{"name":"{name}","stereo":[{}]}}"#, q.join(",")).ok();
                    }
                    _ => { writeln!(w, r#"{{"name":"{name}","error":"unparseable"}}"#).ok(); }
                }
            }
            ExitCode::SUCCESS
        }
        "bondtypes" => {
            let Some(path) = std::env::args().nth(2) else { eprintln!("usage"); return ExitCode::from(2); };
            let f = std::fs::File::open(&path).expect("open corpus");
            let stdout = std::io::stdout(); let mut w = BufWriter::new(stdout.lock());
            for line in std::io::BufReader::new(f).lines() {
                let Ok(line) = line else { break };
                let mut it = line.split_whitespace();
                let (Some(smiles), Some(name)) = (it.next(), it.next()) else { continue };
                match bb_rdkit::bond_types(smiles) {
                    Ok(v) => { let q: Vec<String> = v.iter().map(|x| x.to_string()).collect();
                        writeln!(w, r#"{{"name":"{name}","bt":[{}]}}"#, q.join(",")).ok(); }
                    _ => { writeln!(w, r#"{{"name":"{name}","error":"x"}}"#).ok(); }
                }
            }
            ExitCode::SUCCESS
        }
        "rings" => {
            let Some(path) = std::env::args().nth(2) else {
                eprintln!("usage: bb-uff-dump rings <corpus.smi>");
                return ExitCode::from(2);
            };
            let f = std::fs::File::open(&path).expect("open corpus");
            let stdout = std::io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            for line in std::io::BufReader::new(f).lines() {
                let Ok(line) = line else { break };
                let mut it = line.split_whitespace();
                let (Some(smiles), Some(name)) = (it.next(), it.next()) else { continue };
                match bb_rdkit::sssr_rings(smiles) {
                    Ok(v) => {
                        let q: Vec<String> = v.iter().map(|x| x.to_string()).collect();
                        writeln!(w, r#"{{"name":"{name}","rings":[{}]}}"#, q.join(",")).ok();
                    }
                    _ => { writeln!(w, r#"{{"name":"{name}","error":"unparseable"}}"#).ok(); }
                }
            }
            ExitCode::SUCCESS
        }
        "labels" | "bondlen" | "dist" | "bounds12" | "bounds13" | "bounds14" => {
            let Some(path) = std::env::args().nth(2) else {
                eprintln!("usage: bb-uff-dump {mode} <corpus.smi>");
                return ExitCode::from(2);
            };
            let Ok(f) = std::fs::File::open(&path) else {
                eprintln!("bb-uff-dump: cannot open {path}");
                return ExitCode::from(2);
            };
            let stdout = std::io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            let (mut ok, mut bad) = (0usize, 0usize);
            for line in std::io::BufReader::new(f).lines() {
                let Ok(line) = line else { break };
                let mut it = line.split_whitespace();
                let (Some(smiles), Some(name)) = (it.next(), it.next()) else {
                    continue;
                };
                let numeric = |key: &str, v: Vec<f64>| {
                    let q: Vec<String> = v.iter().map(|x| format!("{x}")).collect();
                    format!(r#""{key}":[{}]"#, q.join(","))
                };
                let row = match mode.as_str() {
                    "labels" => bb_rdkit::uff_atom_labels(smiles).ok().filter(|v| !v.is_empty()).map(|labels| {
                        let q: Vec<String> = labels.iter().map(|l| format!("{l:?}")).collect();
                        format!(r#""labels":[{}]"#, q.join(","))
                    }),
                    "bondlen" => bb_rdkit::uff_bond_rest_lengths(smiles).ok()
                        .filter(|v| !v.is_empty()).map(|v| numeric("bondlen", v)),
                    "dist" => bb_rdkit::topo_distance_matrix(smiles).ok()
                        .filter(|v| !v.is_empty()).map(|v| numeric("dist", v)),
                    // set12 + VDW only (set13=set14=set15=false)
                    "bounds12" => bb_rdkit::raw_bounds_stage(smiles, false, false, false).ok()
                        .filter(|v| !v.is_empty()).map(|v| numeric("bounds12", v)),
                    // set12 + set13 + VDW (set14=set15=false)
                    "bounds13" => bb_rdkit::raw_bounds_stage(smiles, true, false, false).ok()
                        .filter(|v| !v.is_empty()).map(|v| numeric("bounds13", v)),
                    // set12 + set13 + set14 + VDW (set15=false)
                    _ => bb_rdkit::raw_bounds_stage(smiles, true, true, false).ok()
                        .filter(|v| !v.is_empty()).map(|v| numeric("bounds14", v)),
                };
                match row {
                    Some(body) => {
                        writeln!(w, r#"{{"name":"{name}",{body}}}"#).expect("write");
                        ok += 1;
                    }
                    None => {
                        writeln!(w, r#"{{"name":"{name}","error":"unparseable"}}"#).expect("write");
                        bad += 1;
                    }
                }
            }
            w.flush().expect("flush");
            eprintln!("  {ok} molecules ({mode}), {bad} unparseable");
            ExitCode::SUCCESS
        }
        _ => {
            eprintln!("usage: bb-uff-dump <table|labels|bondlen|dist|bounds12|bounds13> [corpus.smi]");
            ExitCode::from(2)
        }
    }
}
