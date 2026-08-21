//! AMSOL7.1 input files.
//!
//! One file per solvent, each a keyword block followed by the molecule name, the atom count, a
//! blank line, and the Z-matrix from [`crate::zmatrix`]:
//!
//! ```text
//! CHARGE=0 AM1 1SCF TLIMIT=15 GEO-OK SM5.42R
//! & SOLVNT=WATER
//! mol0001 42
//!
//! C    0.000000  0 …
//! ```

use crate::config::SolvConfig;
use crate::zmatrix::{render, ZLine};

/// The solvent an input file selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Solvent {
    Water,
    Hexadecane,
}

impl Solvent {
    /// The filename the pipeline writes this solvent's input to.
    pub fn input_filename(self) -> &'static str {
        match self {
            Solvent::Water => "temp.in-wat",
            Solvent::Hexadecane => "temp.in-hex",
        }
    }

    /// The filename the pipeline captures this solvent's AMSOL output to.
    pub fn output_filename(self) -> &'static str {
        match self {
            Solvent::Water => "temp.o-wat",
            Solvent::Hexadecane => "temp.o-hex",
        }
    }

    /// The configured continuation lines for this solvent.
    fn solvent_lines(self, cfg: &SolvConfig) -> &[String] {
        match self {
            Solvent::Water => &cfg.water_solvent,
            Solvent::Hexadecane => &cfg.hexadecane_solvent,
        }
    }
}

/// Build the AMSOL input file for one solvent.
///
/// `charge` is the molecule's net formal charge and `name` its identifier; `zlines` is the
/// Z-matrix, whose length gives the atom count written on the name line.
pub fn amsol_input(
    cfg: &SolvConfig,
    solvent: Solvent,
    name: &str,
    charge: i32,
    zlines: &[ZLine],
) -> String {
    let mut s = String::new();
    s.push_str(&format!("CHARGE={charge} {}\n", cfg.method));
    for line in solvent.solvent_lines(cfg) {
        s.push_str(&format!("& {line}\n"));
    }
    s.push_str(&format!("{name} {}\n", zlines.len()));
    s.push('\n');
    s.push_str(&render(zlines));
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{AmsolFile, SolvConfig};
    use crate::zmatrix::to_zmatrix;

    fn cfg() -> SolvConfig {
        SolvConfig::resolve(&AmsolFile::default(), None).expect("resolved")
    }

    fn methane_zlines() -> Vec<ZLine> {
        // C with four H at tetrahedral positions
        let z = vec![6u8, 1, 1, 1, 1];
        let coords = vec![
            0.0, 0.0, 0.0, 0.63, 0.63, 0.63, -0.63, -0.63, 0.63, -0.63, 0.63, -0.63, 0.63, -0.63,
            -0.63,
        ];
        let bonds = vec![[0, 1], [0, 2], [0, 3], [0, 4]];
        to_zmatrix(&z, &coords, &bonds).expect("zmatrix built")
    }

    #[test]
    fn water_header_matches_the_pipeline_format() {
        let zl = methane_zlines();
        let out = amsol_input(&cfg(), Solvent::Water, "mol0001", 0, &zl);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "CHARGE=0 AM1 1SCF TLIMIT=15 GEO-OK SM5.42R");
        assert_eq!(lines[1], "& SOLVNT=WATER");
        assert_eq!(lines[2], "mol0001 5");
        assert_eq!(lines[3], "");
        assert!(
            lines[4].starts_with("C "),
            "z-matrix follows: {:?}",
            lines[4]
        );
        assert_eq!(lines.len(), 4 + zl.len());
    }

    #[test]
    fn hexadecane_header_carries_both_continuation_lines() {
        let zl = methane_zlines();
        let out = amsol_input(&cfg(), Solvent::Hexadecane, "mol0001", 0, &zl);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "CHARGE=0 AM1 1SCF TLIMIT=15 GEO-OK SM5.42R");
        assert_eq!(
            lines[1],
            "& SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93"
        );
        assert_eq!(lines[2], "& DIELEC=2.06 FACARB=0.00 FEHALO=0.00 DEV");
        assert_eq!(lines[3], "mol0001 5");
        assert_eq!(lines[4], "");
    }

    #[test]
    fn charge_is_written_signed() {
        let zl = methane_zlines();
        for (q, expect) in [(0, "CHARGE=0"), (1, "CHARGE=1"), (-2, "CHARGE=-2")] {
            let out = amsol_input(&cfg(), Solvent::Water, "m", q, &zl);
            assert!(
                out.starts_with(expect),
                "charge {q} rendered as {:?}",
                out.lines().next()
            );
        }
    }

    #[test]
    fn atom_count_tracks_the_zmatrix_length() {
        let zl = methane_zlines();
        let out = amsol_input(&cfg(), Solvent::Water, "m", 0, &zl);
        assert_eq!(out.lines().nth(2), Some("m 5"));
        assert_eq!(zl.len(), 5);
    }

    #[test]
    fn configured_keywords_replace_the_defaults() {
        let explicit = AmsolFile {
            method: Some("AM1 1SCF".into()),
            water_solvent: Some(vec!["SOLVNT=CUSTOM".into(), "EXTRA=1".into()]),
            ..Default::default()
        };
        let c = SolvConfig::resolve(&explicit, None).expect("resolved");
        let out = amsol_input(&c, Solvent::Water, "m", 0, &methane_zlines());
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "CHARGE=0 AM1 1SCF");
        assert_eq!(lines[1], "& SOLVNT=CUSTOM");
        assert_eq!(lines[2], "& EXTRA=1");
    }

    #[test]
    fn filenames_match_the_pipeline() {
        assert_eq!(Solvent::Water.input_filename(), "temp.in-wat");
        assert_eq!(Solvent::Water.output_filename(), "temp.o-wat");
        assert_eq!(Solvent::Hexadecane.input_filename(), "temp.in-hex");
        assert_eq!(Solvent::Hexadecane.output_filename(), "temp.o-hex");
    }
}
