//! Ad-hoc timing harness (dev): times the RDKit-free CPU build stage, broken into sub-stages so it can
//! be compared apples-to-apples with Divya's `build_ligands.py` per-stage report (`rdkit_conf_gen`,
//! `strain`, `db2`). AMSOL solvation is excluded from both sides — it is an identical external
//! subprocess. Pin to one core with `RAYON_NUM_THREADS=1` for a true single-thread comparison.
//!   timeit <input.smi>

use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("usage: timeit <input.smi>");
    let text = std::fs::read_to_string(&path).expect("read input");
    let solv0 = |n: usize| bb_solv::solv::SolvFile {
        name: "fake".into(), charge: 0.0, area: 0.0, tot_diff_pol: 0.0, tot_diff_apol: 0.0,
        tot_diff_pol_plus_apol: 0.0,
        atoms: vec![bb_solv::solv::SolvAtom { charge: 0.0, diff_pol: 0.0, area: 0.0, diff_apol: 0.0, diff_atomic_solv: 0.0 }; n],
    };

    let threads = std::env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "default(all cores)".into());
    eprintln!("RAYON_NUM_THREADS={threads}\n");

    let t0 = Instant::now();
    let (mut ok, mut confs_total) = (0usize, 0usize);
    let (mut te, mut ta, mut ts): (f64, f64, f64) = (0.0, 0.0, 0.0);
    for line in text.lines() {
        let mut it = line.split_whitespace();
        let (Some(smiles), Some(name)) = (it.next(), it.next()) else { continue };

        // conf-gen (perceive + spec + embed 200) — comparable to her `rdkit_conf_gen`
        let e0 = Instant::now();
        let (spec, p) = bb_spec::build_native_perceived(smiles).expect("spec");
        let raw = bb_embed::embed_recipe(&spec, 0x0BADC0DE);
        let confs: Vec<Vec<[f64; 3]>> = raw.iter()
            .map(|c| (0..spec.n_atoms).map(|i| [c.coords[3*i], c.coords[3*i+1], c.coords[3*i+2]]).collect())
            .collect();
        let e_embed = e0.elapsed().as_secs_f64();

        // typing + strain (assemble) — comparable to her SYBYL typing + `strain`
        let a0 = Instant::now();
        let mol = bb_output::assemble::assemble(name, smiles, &p, confs);
        let e_asm = a0.elapsed().as_secs_f64();

        // serialize both members (mol2 + db2) — comparable to her `db2`
        let n = mol.atom_num.len();
        let solv = solv0(n);
        let s0 = Instant::now();
        let _mol2 = bb_output::write_mol2(&mol, &vec![0.0; n], 0);
        let _db2 = bb_db2::write_entry(&bb_output::build(&mol, &solv).expect("db2"));
        let e_ser = s0.elapsed().as_secs_f64();

        te += e_embed; ta += e_asm; ts += e_ser;
        confs_total += mol.atom_xyz.len(); ok += 1;
        eprintln!("{name}: {n} atoms, {} confs | conf_gen={e_embed:.2}s type+strain={e_asm:.2}s serialize={e_ser:.2}s  total={:.2}s",
            mol.atom_xyz.len(), e_embed + e_asm + e_ser);
    }
    let el = t0.elapsed().as_secs_f64();
    eprintln!("\n{ok} molecules, {confs_total} confs [NO AMSOL]");
    eprintln!("  conf_gen(embed):  {te:.2}s  ({:.2}s/mol)", te / ok as f64);
    eprintln!("  type+strain:      {ta:.2}s  ({:.2}s/mol)", ta / ok as f64);
    eprintln!("  serialize mol2+db2: {ts:.2}s  ({:.2}s/mol)", ts / ok as f64);
    eprintln!("  TOTAL CPU stage:  {el:.2}s  ({:.2}s/mol)", el / ok as f64);
}
