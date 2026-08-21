//! Emit bb-bounds' computed quantities per molecule, for parity gates against RDKit.
//!
//!   bb-bounds-dump labels  <corpus.smi>   UFF atom-type labels (post-AddHs atom order)
//!   bb-bounds-dump bondlen <corpus.smi>   UFF bond rest lengths (post-AddHs bond order)
//!   bb-bounds-dump dist    <corpus.smi>   topological distance matrix (flat n*n row-major)
//!
//! One JSON object per line, keyed by corpus name, matching the shapes `bb-uff-dump` captures from
//! RDKit so the gate compares directly. A SMILES bb-perceive cannot parse is emitted with an
//! "error" field rather than dropped, so the count matches the corpus.

use std::io::{BufRead, BufWriter, Write};
use std::process::ExitCode;

fn main() -> ExitCode {
    let mode = std::env::args().nth(1).unwrap_or_default();
    if !matches!(
        mode.as_str(),
        "labels"
            | "bondlen"
            | "dist"
            | "bounds12"
            | "bounds13"
            | "bounds14"
            | "bounds15"
            | "rings"
            | "arom"
            | "stereo"
            | "chiraltags"
    ) {
        eprintln!("usage: bb-bounds-dump <labels|bondlen|dist|bounds12|bounds13|bounds14|bounds15|stereo|chiraltags> <corpus.smi>");
        return ExitCode::from(2);
    }
    let Some(path) = std::env::args().nth(2) else {
        eprintln!("usage: bb-bounds-dump {mode} <corpus.smi>");
        return ExitCode::from(2);
    };
    let Ok(f) = std::fs::File::open(&path) else {
        eprintln!("bb-bounds-dump: cannot open {path}");
        return ExitCode::from(2);
    };

    let out = std::io::stdout();
    let mut w = BufWriter::new(out.lock());
    for line in std::io::BufReader::new(f).lines() {
        let Ok(line) = line else { break };
        let mut it = line.split_whitespace();
        let (Some(smiles), Some(name)) = (it.next(), it.next()) else {
            continue;
        };
        match bb_perceive::smarts_match::perceive(smiles) {
            Ok(mol) => {
                let body = match mode.as_str() {
                    "labels" => {
                        let q: Vec<String> = bb_bounds::uff::atom_labels(&mol)
                            .iter()
                            .map(|l| format!("{l:?}"))
                            .collect();
                        format!(r#""labels":[{}]"#, q.join(","))
                    }
                    "bondlen" => {
                        let q: Vec<String> = bb_bounds::uff::bond_rest_lengths(&mol)
                            .iter()
                            .map(|x| format!("{x}"))
                            .collect();
                        format!(r#""bondlen":[{}]"#, q.join(","))
                    }
                    "dist" => {
                        let d =
                            bb_bounds::dist::distance_matrix(mol.atomic_numbers.len(), &mol.bonds);
                        let q: Vec<String> = d.iter().map(|x| format!("{x}")).collect();
                        format!(r#""dist":[{}]"#, q.join(","))
                    }
                    "bounds12" => {
                        let b = bb_bounds::bounds::bounds_set12_vdw(&mol);
                        let q: Vec<String> = b.iter().map(|x| format!("{x}")).collect();
                        format!(r#""bounds12":[{}]"#, q.join(","))
                    }
                    "bounds13" => {
                        let b = bb_bounds::bounds::topol_bounds(&mol, true, false, false);
                        let q: Vec<String> = b.iter().map(|x| format!("{x}")).collect();
                        format!(r#""bounds13":[{}]"#, q.join(","))
                    }
                    "bounds14" => {
                        let b = bb_bounds::bounds::topol_bounds(&mol, true, true, false);
                        let q: Vec<String> = b.iter().map(|x| format!("{x}")).collect();
                        format!(r#""bounds14":[{}]"#, q.join(","))
                    }
                    "bounds15" => {
                        // the full setTopolBounds, compared against the raw_bounds fixture (key "bounds")
                        let b = bb_bounds::bounds::bounds_full(&mol);
                        let q: Vec<String> = b.iter().map(|x| format!("{x}")).collect();
                        format!(r#""bounds":[{}]"#, q.join(","))
                    }
                    "arom" => {
                        // perceived per-atom then per-bond aromaticity, matching bb-uff-dump's
                        // aromatic_perception oracle: [n_atoms, atomArom.., (begin,end,bondArom)..]
                        let mut q: Vec<String> = vec![mol.atomic_numbers.len().to_string()];
                        q.extend(mol.aromatic.iter().map(|a| i32::from(*a).to_string()));
                        for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
                            q.push(a.to_string());
                            q.push(b.to_string());
                            q.push(i32::from(mol.bond_aromatic[bi]).to_string());
                        }
                        format!(r#""arom":[{}]"#, q.join(","))
                    }
                    "rings" => {
                        // symmetrized SSSR, flat with -1 ring separators, matching bb-uff-dump rings
                        let mut q: Vec<String> = Vec::new();
                        for r in &mol.rings {
                            q.extend(r.iter().map(|a| a.to_string()));
                            q.push("-1".to_string());
                        }
                        format!(r#""rings":[{}]"#, q.join(","))
                    }
                    "chiraltags" => {
                        use bb_perceive::chirality::ChiralTag;
                        // ints matching RDKit ChiralType: UNSPECIFIED=0, CW=1, CCW=2
                        let q: Vec<String> = mol
                            .chiral_tags
                            .iter()
                            .map(|t| {
                                match t {
                                    ChiralTag::None => "0",
                                    ChiralTag::Cw => "1",
                                    ChiralTag::Ccw => "2",
                                }
                                .to_string()
                            })
                            .collect();
                        format!(r#""tags":[{}]"#, q.join(","))
                    }
                    _ => {
                        use bb_perceive::smarts_match::Stereo;
                        // enum ints matching RDKit: NONE=0 ANY=1 Z=2 E=3 CIS=4 TRANS=5
                        let q: Vec<String> = mol
                            .bond_stereo
                            .iter()
                            .zip(&mol.stereo_atoms)
                            .flat_map(|(s, a)| {
                                let e = match s {
                                    Stereo::None => 0,
                                    Stereo::Any => 1,
                                    Stereo::Z => 2,
                                    Stereo::E => 3,
                                    Stereo::Cis => 4,
                                    Stereo::Trans => 5,
                                };
                                [e.to_string(), a[0].to_string(), a[1].to_string()]
                            })
                            .collect();
                        format!(r#""stereo":[{}]"#, q.join(","))
                    }
                };
                writeln!(w, r#"{{"name":"{name}",{body}}}"#).ok();
            }
            Err(e) => {
                writeln!(w, r#"{{"name":"{name}","error":"{e}"}}"#).ok();
            }
        }
    }
    ExitCode::SUCCESS
}
