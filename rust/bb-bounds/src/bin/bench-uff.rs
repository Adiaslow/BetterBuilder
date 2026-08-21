//! Time UFF atom typing with the same discipline as the C++ baseline.
//!
//!   bench-uff <corpus.smi>
//!
//! Mirrors `validation/bench/time_cpp_uff.cc`: molecules are perceived up front (outside the timer),
//! then only the typing component — `atom_label` + `params_for_label` per atom — is timed, after a
//! warm-up. The equal-work figures (atoms_typed, checksum = sum of r1) must match the C++ bench, or
//! the two are not doing the same work.

use std::time::Instant;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: bench-uff <corpus.smi>");
    let text = std::fs::read_to_string(&path).expect("read corpus");

    // perceive every molecule up front, outside the timer (the C++ side prepares mols up front too)
    let perceived: Vec<_> = text
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter_map(|smi| bb_perceive::smarts_match::perceive(smi).ok())
        .collect();

    let type_one =
        |mol: &bb_perceive::smarts_match::Perceived, atoms: &mut i64, checksum: &mut f64| {
            for label in bb_bounds::uff::atom_labels(mol) {
                if let Some(p) = bb_bounds::uff::params_for_label(&label) {
                    *checksum += p.r1;
                    *atoms += 1;
                }
            }
        };

    for mol in perceived.iter().take(5) {
        let (mut a, mut s) = (0i64, 0.0f64);
        type_one(mol, &mut a, &mut s);
    }

    let mut atoms = 0i64;
    let mut checksum = 0.0f64;
    let t0 = Instant::now();
    for mol in &perceived {
        type_one(mol, &mut atoms, &mut checksum);
    }
    let secs = t0.elapsed().as_secs_f64();

    eprintln!(
        "rust bb-bounds UFF typing (perceived mols): {} mols in {secs:.5} s = {:.5} ms/mol, atoms_typed {atoms}, checksum {checksum:.12}",
        perceived.len(),
        1000.0 * secs / perceived.len() as f64
    );
}
