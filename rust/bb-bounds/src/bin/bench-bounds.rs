//! Time the full setTopolBounds with the same discipline as the C++ baseline.
//!
//!   bench-bounds <corpus.smi>
//!
//! Mirrors `validation/bench/time_cpp_bounds.cc`: molecules are perceived up front (outside the
//! timer), then only the bounds-matrix construction (`bounds_full`) is timed, after a warm-up. The
//! equal-work figures — `constrained_pairs` (upper bounds tightened below MAX_UPPER, an exact
//! integer) and `checksum` (sum of all cells) — must match the C++ bench, or the two are not doing
//! the same work.

use std::time::Instant;

const MAX_UPPER: f64 = 1000.0;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: bench-bounds <corpus.smi>");
    let text = std::fs::read_to_string(&path).expect("read corpus");
    let perceived: Vec<_> = text
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter_map(|smi| bb_perceive::smarts_match::perceive(smi).ok())
        .collect();

    let build_one =
        |mol: &bb_perceive::smarts_match::Perceived, constrained: &mut i64, checksum: &mut f64| {
            let b = bb_bounds::bounds::bounds_full(mol);
            let n = mol.atomic_numbers.len();
            for i in 0..n {
                for j in (i + 1)..n {
                    if b[i * n + j] < MAX_UPPER {
                        *constrained += 1;
                    }
                }
            }
            for &x in &b {
                *checksum += x;
            }
        };

    for mol in perceived.iter().take(5) {
        let (mut c, mut s) = (0i64, 0.0f64);
        build_one(mol, &mut c, &mut s);
    }

    let mut constrained = 0i64;
    let mut checksum = 0.0f64;
    let t0 = Instant::now();
    for mol in &perceived {
        build_one(mol, &mut constrained, &mut checksum);
    }
    let secs = t0.elapsed().as_secs_f64();

    eprintln!(
        "rust bb-bounds setTopolBounds (perceived mols): {} mols in {secs:.5} s = {:.5} ms/mol, constrained_pairs {constrained}, checksum {checksum:.6}",
        perceived.len(),
        1000.0 * secs / perceived.len() as f64
    );
}
