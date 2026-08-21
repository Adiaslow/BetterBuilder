//! bb-solv — solvation input preparation for AMSOL7.1.
//!
//! [`zmatrix`] converts conformer coordinates to MOPAC internal coordinates, [`input`] wraps those
//! in the per-solvent AMSOL input files, [`config`] resolves the AMSOL executable and the keyword
//! blocks those files carry, [`run`] executes AMSOL over them, and [`solv`] turns its output into
//! the `.solv` file the DOCK pipeline consumes.

pub mod config;
pub mod input;
pub mod run;
pub mod solv;
pub mod zmatrix;
