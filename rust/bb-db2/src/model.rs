//! The db2 conformer-hierarchy data model — one [`Db2Entry`] per stacked molecule entry.
//!
//! A db2 file is one or more entries, each an atom **pool** (`A`/`B`) plus a coordinate **hierarchy**:
//! a rigid core (`R`), every flexible-atom position (`X`), conformations grouping coordinate runs
//! (`C`), and sets selecting conformations per input conformer (`S`). The field set and record order
//! are pinned by `dock3/src/db2formats.f` (the reader) and the `mol2db2` writer; see
//! `docs/…/09-db2-format.md`. Numeric ids are 1-based on disk and kept 1-based here so a parsed entry
//! re-serializes byte-for-byte.

/// One `A` record: an atom in the pool (invariant across conformers).
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    pub num: i32,
    pub name: String,
    pub sybyl: String,
    pub vdwtype: i32,
    pub color: i32,
    pub charge: f64,
    pub polar_solv: f64,
    pub apolar_solv: f64,
    pub total_solv: f64,
    pub surface: f64,
}

/// One `B` record: a bond between two 1-based pool atoms; `btype` is a SYBYL-ish code (`1`, `2`, `ar`).
#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    pub num: i32,
    pub a1: i32,
    pub a2: i32,
    pub btype: String,
}

/// One `X` record: a single atom's position within a single conformation.
#[derive(Clone, Debug, PartialEq)]
pub struct Coord {
    pub num: i32,
    pub atom: i32,
    pub conf: i32,
    pub xyz: [f64; 3],
}

/// One `R` record: a rigid-core atom's single position (with its matching-sphere color).
#[derive(Clone, Debug, PartialEq)]
pub struct Rigid {
    pub num: i32,
    pub color: i32,
    pub xyz: [f64; 3],
}

/// One `C` record: conformation `num` is the inclusive 1-based `X` range `[start, end]`.
#[derive(Clone, Debug, PartialEq)]
pub struct Conf {
    pub num: i32,
    pub start: i32,
    pub end: i32,
}

/// One `S` set: an input conformer as a header plus member lines of conformation indices.
///
/// `broken` = clash-check failure flag, `hydrogens` = terminal-H treatment (0 input / 1 reset /
/// 2 rotated), and the two strain values come from the strain stage. `member_lines` preserves the
/// on-disk line grouping (up to 8 indices per line) so the `#lines` count and layout round-trip.
#[derive(Clone, Debug, PartialEq)]
pub struct Set {
    pub num: i32,
    pub broken: i32,
    pub hydrogens: i32,
    pub total_strain: f64,
    pub max_strain: f64,
    pub member_lines: Vec<Vec<i32>>,
}

/// One complete molecule entry (`M`…`E`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Db2Entry {
    /// `T` color-table entries: `(code, name)`, only the non-default colors.
    pub color_table: Vec<(i32, String)>,
    // M1 identity
    pub name: String,
    pub protname: String,
    // M2 solvation totals
    pub total_charge: f64,
    pub total_polar_solv: f64,
    pub total_apolar_solv: f64,
    pub total_solv: f64,
    pub surface_area: f64,
    // M3 / M4 / M5
    pub smiles: String,
    pub longname: String,
    pub m5: f64,
    /// Any `M` lines beyond the canonical five (rare), stored raw (without the leading `M`).
    pub extra_m: Vec<String>,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub coords: Vec<Coord>,
    pub rigid: Vec<Rigid>,
    pub confs: Vec<Conf>,
    pub sets: Vec<Set>,
    /// `D` cluster block lines (rare — most builds emit none), stored raw for faithful round-trip.
    pub cluster_lines: Vec<String>,
}
