//! Self-contained format guards. The full byte-exactness proof is `bb-db2-roundtrip` over the real
//! DOCK-pipeline db2 corpus (3287 files, all byte-identical); these tests lock the format in-repo
//! without that external data: an exact-bytes assertion per record type (hand-verified against real
//! lines) and a write→parse→write idempotence check.

use bb_db2::model::*;
use bb_db2::{parse, write_file};

/// A small entry exercising every record type. Counts in M1 / S headers are recomputed by the
/// writer, so they always agree with the records here.
fn sample() -> Db2Entry {
    Db2Entry {
        color_table: vec![],
        name: "LF2B76".into(),
        protname: "fake".into(),
        total_charge: 0.0,
        total_polar_solv: -13.63,
        total_apolar_solv: 6.1,
        total_solv: -7.53,
        surface_area: 441.63,
        smiles: "COCC1CCN(C(=O)C2(C)CN(C(=O)C3CCC(C(F)(F)F)OC3)CCO2)CC1".into(),
        longname: "fake".into(),
        m5: 999.999,
        extra_m: vec![],
        atoms: vec![Atom {
            num: 1,
            name: "C1".into(),
            sybyl: "C.3".into(),
            vdwtype: 5,
            color: 7,
            charge: 0.03,
            polar_solv: -0.03,
            apolar_solv: -0.02,
            total_solv: -0.05,
            surface: 11.01,
        }],
        bonds: vec![Bond { num: 1, a1: 1, a2: 2, btype: "1".into() }],
        coords: vec![Coord { num: 1, atom: 3, conf: 1, xyz: [-1.4054, -2.5097, -0.3192] }],
        rigid: vec![Rigid { num: 1, color: 7, xyz: [-1.4054, -2.5097, -0.3192] }],
        confs: vec![Conf { num: 1, start: 1, end: 17 }],
        sets: vec![Set {
            num: 1,
            broken: 1,
            hydrogens: 0,
            total_strain: 3.297,
            max_strain: 1.602,
            member_lines: vec![vec![1, 2, 3, 5, 30, 31, 32, 33], vec![11, 18, 19]],
        }],
        cluster_lines: vec![],
    }
}

#[test]
fn exact_bytes_per_record() {
    let out = write_file(&[sample()]);
    let lines: Vec<&str> = out.lines().collect();
    // hand-verified against real db2 lines (LF2B76.0.N.db2)
    assert_eq!(lines[0], "M           LF2B76      fake   1   1      1      1      1      1      5      0");
    assert_eq!(lines[1], "M   +0.0000    -13.630     +6.100     -7.530   441.630");
    assert_eq!(lines[2], "M COCC1CCN(C(=O)C2(C)CN(C(=O)C3CCC(C(F)(F)F)OC3)CCO2)CC1                      ");
    assert_eq!(lines[4], "M  +999.9990");
    assert_eq!(lines[5], "A   1 C1   C.3    5  7   +0.0300     -0.030     -0.020     -0.050    11.010");
    assert_eq!(lines[6], "B   1   1   2 1 ");
    assert_eq!(lines[7], "X         1   3      1   -1.4054   -2.5097   -0.3192");
    assert_eq!(lines[8], "R      1  7   -1.4054   -2.5097   -0.3192");
    assert_eq!(lines[9], "C      1         1        17");
    assert_eq!(lines[10], "S      1      2  11 1 0      +3.297      +1.602");
    assert_eq!(lines[11], "S      1      1 8      1      2      3      5     30     31     32     33");
    assert_eq!(lines[12], "S      1      2 3     11     18     19");
    assert_eq!(lines[13], "E");
}

#[test]
fn write_parse_write_is_idempotent() {
    let first = write_file(&[sample()]);
    let parsed = parse(&first).expect("parse");
    let second = write_file(&parsed);
    assert_eq!(first, second, "write→parse→write must be a fixed point");
}

#[test]
fn multi_entry_stacks() {
    let stacked = write_file(&[sample(), sample()]);
    assert_eq!(stacked.lines().filter(|l| *l == "E").count(), 2);
    let parsed = parse(&stacked).expect("parse");
    assert_eq!(parsed.len(), 2);
    assert_eq!(write_file(&parsed), stacked);
}
