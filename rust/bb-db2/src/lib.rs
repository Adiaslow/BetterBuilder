//! bb-db2 — the DOCK `.db2` conformer-hierarchy format: a data model ([`model`]), a faithful
//! [`parse`]r, and a byte-exact [`write`]r.
//!
//! The `.db2` file is the last artifact of the ligand-build pipeline — the atom pool plus the
//! rigid/flexible coordinate hierarchy DOCK docks against. The format is pinned by
//! `dock3/src/db2formats.f` (the canonical Fortran read side) and the `mol2db2` writer's printf
//! templates. This crate reproduces the *format* exactly (proven by round-tripping real
//! pipeline-built db2 files to byte identity); assembling an entry from BetterBuilder's own
//! conformers + solvation is a separate, later stage.

pub mod model;
pub mod parse;
pub mod write;

pub use model::Db2Entry;
pub use parse::{parse, ParseError};
pub use write::{write_entry, write_file};
