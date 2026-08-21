//! Emit assigned torsions per molecule, for comparison against RDKit.
//!
//!     bb-torsions <corpus.smi>
use std::io::{BufRead, BufWriter, Write};

fn main() {
    let path = std::env::args().nth(1).expect("usage: bb-torsions <corpus.smi>");
    let lib = bb_perceive::torsion_lib::patterns();
    let asts = bb_perceive::torsions::compile(&lib);
    let f = std::fs::File::open(&path).expect("open corpus");
    let out = std::io::stdout();
    let mut w = BufWriter::new(out.lock());
    for line in std::io::BufReader::new(f).lines() {
        let line = line.expect("read");
        let mut it = line.split_whitespace();
        let (Some(smiles), Some(name)) = (it.next(), it.next()) else { continue };
        match bb_perceive::smarts_match::perceive(smiles) {
            Ok(mol) => {
                let (mut ts, mut done) = bb_perceive::torsions::assign(&mol, &lib, &asts);
                let (bk, _imp) = bb_perceive::torsions::basic_knowledge(&mol, &mut done);
                ts.extend(bk);
                let items: Vec<String> = ts
                    .iter()
                    .map(|t| {
                        let v: Vec<String> = t.v.iter().map(|x| format!("{x:.6}")).collect();
                        format!(
                            r#"{{"atoms":[{},{},{},{}],"v":[{}]}}"#,
                            t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3], v.join(",")
                        )
                    })
                    .collect();
                writeln!(w, r#"{{"name":"{name}","torsions":[{}]}}"#, items.join(",")).ok();
            }
            Err(e) => {
                writeln!(w, r#"{{"name":"{name}","error":"{e}"}}"#).ok();
            }
        }
    }
}
