//! Byte-exact db2 serialization: [`Db2Entry`] → the `T/M/A/B/X/R/C/S/D/E` text a DOCK build reads.
//!
//! Each record reproduces the `mol2db2` writer's printf template verbatim (documented in
//! `09-db2-format.md`, consistent with `dock3/src/db2formats.f`). The template *includes* the tag
//! letter and the literal separating spaces, so a Rust `format!` mirrors it field-for-field:
//! `%3d`→`{:3}`, `%-4s`→`{:<4}`, `%+9.4f`→`{:+9.4}`, `%9.3f`→`{:9.3}`. Counts and per-set line
//! tallies are recomputed from the data so they always agree with the records that follow.

use crate::model::Db2Entry;
use std::fmt::Write as _;

/// The last `n` bytes of `s` (db2 strings are ASCII), reproducing Python's `s[-n:]` truncation the
/// `mol2db2` writer applies before field-formatting (`name[-16:]`, `smiles[-76:]`, …).
fn last(s: &str, n: usize) -> &str {
    let len = s.len();
    if len > n {
        &s[len - n..]
    } else {
        s
    }
}

/// Serialize one entry (no trailing newline beyond the final `E` line's).
pub fn write_entry(e: &Db2Entry) -> String {
    let mut s = String::new();
    write_entry_into(&mut s, e);
    s
}

/// Serialize a stack of entries into one db2 file body.
pub fn write_file(entries: &[Db2Entry]) -> String {
    let mut s = String::new();
    for e in entries {
        write_entry_into(&mut s, e);
    }
    s
}

fn write_entry_into(s: &mut String, e: &Db2Entry) {
    // T — color table (only non-default colors; usually none); colorName is right-justified
    for (code, name) in &e.color_table {
        let _ = writeln!(s, "T {code:2} {name:>8}");
    }

    // M1 — identity + counts (counts recomputed from the records below)
    // M1: name/protname right-justified, truncated to their last 16/9 chars; #maxmlines is a literal 5.
    let _ = writeln!(
        s,
        "M {:>16} {:>9} {:3} {:3} {:6} {:6} {:6} {:6} {:6} {:6}",
        last(&e.name, 16),
        last(&e.protname, 9),
        e.atoms.len(),
        e.bonds.len(),
        e.coords.len(),
        e.confs.len(),
        e.sets.len(),
        e.rigid.len(),
        5,
        e.cluster_lines_count(),
    );
    // M2 — solvation totals
    let _ = writeln!(
        s,
        "M {:+9.4} {:+10.3} {:+10.3} {:+10.3} {:9.3}",
        e.total_charge, e.total_polar_solv, e.total_apolar_solv, e.total_solv, e.surface_area
    );
    // M3 — SMILES, M4 — long name (left-justified in 76, truncated to last 76 chars)
    let _ = writeln!(s, "M {:<76}", last(&e.smiles, 76));
    let _ = writeln!(s, "M {:<76}", last(&e.longname, 76));
    // M5 — arbitrary placeholder
    let _ = writeln!(s, "M {:+10.4}", e.m5);
    // any extra M lines beyond the canonical five, verbatim
    for m in &e.extra_m {
        let _ = writeln!(s, "M {m}");
    }

    // A — atoms
    for a in &e.atoms {
        let _ = writeln!(
            s,
            "A {:3} {:<4} {:<5} {:2} {:2} {:+9.4} {:+10.3} {:+10.3} {:+10.3} {:9.3}",
            a.num, a.name, a.sybyl, a.vdwtype, a.color, a.charge, a.polar_solv, a.apolar_solv,
            a.total_solv, a.surface
        );
    }
    // B — bonds
    for b in &e.bonds {
        let _ = writeln!(s, "B {:3} {:3} {:3} {:<2}", b.num, b.a1, b.a2, b.btype);
    }
    // X — flexible coordinates
    for c in &e.coords {
        let _ = writeln!(
            s,
            "X {:9} {:3} {:6} {:+9.4} {:+9.4} {:+9.4}",
            c.num, c.atom, c.conf, c.xyz[0], c.xyz[1], c.xyz[2]
        );
    }
    // R — rigid-core atoms
    for r in &e.rigid {
        let _ = writeln!(
            s,
            "R {:6} {:2} {:+9.4} {:+9.4} {:+9.4}",
            r.num, r.color, r.xyz[0], r.xyz[1], r.xyz[2]
        );
    }
    // C — conformations
    for c in &e.confs {
        let _ = writeln!(s, "C {:6} {:9} {:9}", c.num, c.start, c.end);
    }
    // S — sets: header then member lines
    for set in &e.sets {
        let nlines = set.member_lines.len();
        let nconfs: usize = set.member_lines.iter().map(Vec::len).sum();
        let _ = writeln!(
            s,
            "S {:6} {:6} {:3} {:1} {:1} {:+11.3} {:+11.3}",
            set.num, nlines, nconfs, set.broken, set.hydrogens, set.total_strain, set.max_strain
        );
        for (li, line) in set.member_lines.iter().enumerate() {
            let _ = write!(s, "S {:6} {:6} {:1}", set.num, li + 1, line.len());
            for conf in line {
                let _ = write!(s, " {conf:6}");
            }
            s.push('\n');
        }
    }
    // D — cluster block, verbatim (rare)
    for line in &e.cluster_lines {
        let _ = writeln!(s, "{line}");
    }
    // E — end of entry
    s.push_str("E\n");
}

impl Db2Entry {
    fn cluster_lines_count(&self) -> usize {
        // #clusters in M1 is the count of `D` header records; we reproduce the raw block, and every
        // build in scope emits none. Count header lines (those beginning with 'D') if present.
        self.cluster_lines
            .iter()
            .filter(|l| l.starts_with('D'))
            .count()
    }
}
