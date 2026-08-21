//! bb-db2-roundtrip — prove the db2 writer is byte-exact against real pipeline output.
//!
//!   bb-db2-roundtrip <file.db2 | dir> [...]
//!
//! For each `.db2` file: read → [`bb_db2::parse`] → [`bb_db2::write_file`] → compare bytes to the
//! original. Real DOCK-pipeline db2 files are the oracle for the *format*; a byte-identical
//! round-trip certifies the serializer reproduces `db2formats.f` / the `mol2db2` writer exactly.
//! Prints a per-file PASS/FAIL and, on the first failure, the first differing line. Exit 0 iff every
//! file round-trips.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

fn collect(path: &Path, out: &mut Vec<PathBuf>) {
    if path.is_dir() {
        if let Ok(rd) = std::fs::read_dir(path) {
            let mut kids: Vec<PathBuf> = rd.flatten().map(|e| e.path()).collect();
            kids.sort();
            for k in kids {
                collect(&k, out);
            }
        }
    } else if path.extension().is_some_and(|x| x == "db2") {
        out.push(path.to_path_buf());
    }
}

/// Return the 1-based number and text of the first differing line, or None if identical.
fn first_diff<'a>(a: &'a str, b: &'a str) -> Option<(usize, &'a str, &'a str)> {
    let mut al = a.lines();
    let mut bl = b.lines();
    let mut n = 0;
    loop {
        n += 1;
        match (al.next(), bl.next()) {
            (None, None) => return None,
            (x, y) if x != y => return Some((n, x.unwrap_or("<EOF>"), y.unwrap_or("<EOF>"))),
            _ => {}
        }
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: bb-db2-roundtrip <file.db2 | dir> [...]");
        return ExitCode::from(2);
    }
    let mut files = Vec::new();
    for a in &args {
        collect(Path::new(a), &mut files);
    }

    let (mut ok, mut fail, mut err) = (0usize, 0usize, 0usize);
    let mut shown = 0;
    for f in &files {
        let orig = match std::fs::read_to_string(f) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("  READ-ERR {}: {e}", f.display());
                err += 1;
                continue;
            }
        };
        let entries = match bb_db2::parse(&orig) {
            Ok(e) => e,
            Err(e) => {
                eprintln!("  PARSE-ERR {}: {e}", f.display());
                err += 1;
                continue;
            }
        };
        let round = bb_db2::write_file(&entries);
        if round == orig {
            ok += 1;
        } else {
            fail += 1;
            if shown < 5 {
                shown += 1;
                let n_orig = orig.lines().count();
                let n_round = round.lines().count();
                eprintln!(
                    "  FAIL {} ({} entries; {n_orig} orig lines, {n_round} round lines)",
                    f.display(),
                    entries.len()
                );
                if let Some((n, o, r)) = first_diff(&orig, &round) {
                    eprintln!("    line {n}:");
                    eprintln!("      orig : {o:?}");
                    eprintln!("      round: {r:?}");
                }
            }
        }
    }
    println!(
        "db2 round-trip: {ok} byte-identical, {fail} differ, {err} unreadable ({} files)",
        files.len()
    );
    if fail == 0 && err == 0 && ok > 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
