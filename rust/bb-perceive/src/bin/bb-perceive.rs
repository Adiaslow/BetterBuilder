//! Emit the heavy-atom graph yowl reads from each SMILES in a corpus, as JSON lines.
//!
//!     bb-perceive <corpus.smi>
use std::io::{BufRead, BufWriter, Write};

fn main() {
    let path = std::env::args().nth(1).expect("usage: bb-perceive <corpus.smi>");
    let f = std::fs::File::open(&path).expect("open corpus");
    let out = std::io::stdout();
    let mut w = BufWriter::new(out.lock());
    for line in std::io::BufReader::new(f).lines() {
        let line = line.expect("read line");
        let mut it = line.split_whitespace();
        let (Some(smiles), Some(name)) = (it.next(), it.next()) else {
            continue;
        };
        match bb_perceive::parse(smiles) {
            Ok(g) => {
                let bonds: Vec<String> =
                    g.bonds.iter().map(|(a, b)| format!("[{a},{b}]")).collect();
                let hc = bb_perceive::hydrogens::counts(&g);
                let ex = bb_perceive::addhs::add_hs(&g, &hc);
                let (conj, hyb) = bb_perceive::hybrid::perceive(&g, &hc);
                // heavy atoms carry hybridization; appended hydrogens stay unspecified
                let mut hyb_s: Vec<String> =
                    hyb.iter().map(|h| (*h as i32).to_string()).collect();
                hyb_s.resize(ex.atomic_numbers.len(), "0".to_string());
                let mut conj_s: Vec<String> =
                    conj.iter().map(|c| i32::from(*c).to_string()).collect();
                conj_s.resize(ex.bonds.len(), "0".to_string());
                // aromaticity as the SMILES writes it, extended over the appended hydrogens
                let mut arom: Vec<i32> =
                    g.atoms.iter().map(|a| i32::from(a.aromatic_as_written)).collect();
                arom.resize(ex.atomic_numbers.len(), 0);
                let exz: Vec<String> = ex.atomic_numbers.iter().map(u8::to_string).collect();
                let exb: Vec<String> =
                    ex.bonds.iter().map(|(a, b)| format!("[{a},{b}]")).collect();
                let js = |v: &[i32]| {
                    v.iter().map(i32::to_string).collect::<Vec<_>>().join(",")
                };
                let rings = bb_perceive::sssr::symmetrized_sssr(g.atoms.len(), &g.bonds);
                let exo = bb_perceive::recipe::count_exo_rotatable(&g, &rings);
                let (seeds, side) = bb_perceive::recipe::recipe(exo);
                let core = bb_perceive::recipe::core_atoms(&rings);
                // integer bond orders over the post-AddHs bond list; appended H bonds are single
                let mut border: Vec<u8> = g
                    .bond_kinds
                    .iter()
                    .enumerate()
                    .map(|(bi, k)| {
                        let (x, y) = g.bonds[bi];
                        let ar = g.atoms[x].aromatic_as_written && g.atoms[y].aromatic_as_written;
                        let c = bb_perceive::hydrogens::contribution(*k, ar);
                        if (c - 2.0).abs() < 1e-9 { 2 } else if (c - 3.0).abs() < 1e-9 { 3 } else { 1 }
                    })
                    .collect();
                border.resize(ex.bonds.len(), 1);
                let mut degs = vec![0usize; ex.atomic_numbers.len()];
                for &(x, y) in &ex.bonds { degs[x] += 1; degs[y] += 1; }
                let angs = bb_perceive::angles::collect(&ex, &border, &degs);
                let angs_s: Vec<String> = angs
                    .iter()
                    .map(|a| format!("[{},{},{},{}]", a.atoms[0], a.atoms[1], a.atoms[2],
                                     i32::from(a.triple)))
                    .collect();
                let core_s: Vec<String> = core.iter().map(usize::to_string).collect();
                let rings_s: Vec<String> = rings
                    .iter()
                    .map(|r| {
                        let a: Vec<String> = r.iter().map(usize::to_string).collect();
                        format!("[{}]", a.join(","))
                    })
                    .collect();
                let z: Vec<String> =
                    g.atoms.iter().map(|a| a.atomic_number.to_string()).collect();
                writeln!(
                    w,
                    r#"{{"name":"{}","n_heavy":{},"z":[{}],"bonds":[{}],"rings":[{}],"explicit_valence":[{}],"implicit_valence":[{}],"total_num_hs":[{}],"addhs_z":[{}],"addhs_bonds":[{}],"arom_written":[{}],"exo":{},"core_seeds":{},"sidechain_confs":{},"core_atoms":[{}],"angles":[{}],"hybridization":[{}],"conjugated":[{}]}}"#,
                    name,
                    g.atoms.len(),
                    z.join(","),
                    bonds.join(","),
                    rings_s.join(","),
                    js(&hc.explicit_valence),
                    js(&hc.implicit_valence),
                    js(&hc.total_num_hs),
                    exz.join(","),
                    exb.join(","),
                    js(&arom),
                    exo,
                    seeds,
                    side,
                    core_s.join(","),
                    angs_s.join(","),
                    hyb_s.join(","),
                    conj_s.join(",")
                )
                .ok();
            }
            Err(e) => {
                writeln!(w, r#"{{"name":"{name}","error":"{e}"}}"#).ok();
            }
        }
    }
}
