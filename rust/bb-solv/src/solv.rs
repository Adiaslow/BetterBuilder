//! AMSOL output → the `.solv` file the DOCK pipeline consumes.
//!
//! [`parse`] reads one AMSOL7.1 output file; [`combine`] takes the water and hexadecane results for
//! the same molecule and forms the water-minus-hexadecane differences; [`SolvFile::render`] writes
//! them in the `.solv` layout.
//!
//! Header: `name n_atoms charge  tot_diff_pol  area  tot_diff_apol  tot_diff_pol_plus_apol`, then
//! one row per atom: `charge  diff_pol  area  diff_apol  diff_atomic_solv`. Charge and area come
//! from the hexadecane run; the energy columns are differences.

use std::fmt;

/// One row of the per-atom table in an AMSOL output file.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomRow {
    pub index: usize,
    pub symbol: String,
    /// CM2 partial charge.
    pub charge: f64,
    /// Atomic polar contribution to the solvation free energy, kcal/mol.
    pub g_p: f64,
    /// Solvent-accessible area, Å².
    pub area: f64,
    /// Surface tension coefficient, kcal/Å².
    pub sigma: f64,
    /// Atomic apolar contribution to the solvation free energy, kcal/mol.
    pub g_cds: f64,
    /// Polar plus apolar for this atom, kcal/mol.
    pub subtotal: f64,
    pub m_value: f64,
}

/// The `Total:` line's five numeric fields.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Totals {
    pub charge: f64,
    pub g_p: f64,
    pub area: f64,
    pub g_cds: f64,
    pub subtotal: f64,
}

/// One parsed AMSOL output file.
#[derive(Clone, Debug, PartialEq)]
pub struct AmsolOutput {
    pub name: String,
    /// The atom count AMSOL echoes from its input.
    pub n_atoms: usize,
    pub atoms: Vec<AtomRow>,
    pub totals: Totals,
}

/// The line AMSOL prints when it gives up. It exits 0 even then, so the output text is the only
/// signal that the run failed.
const AMSOL_FAILURE_MARKER: &str = "was not completed successfully";

#[derive(Debug, PartialEq)]
pub enum SolvError {
    /// AMSOL reported that it could not complete the calculation, carrying its own diagnostic —
    /// commonly a compiled-in size limit such as `MAXLIT` for hydrogen count.
    AmsolFailed { message: String },
    /// No `<name> <count>` line was found.
    NoNameLine,
    /// No `Total:` line was found.
    NoTotalLine,
    /// The per-atom table length disagrees with the echoed atom count.
    AtomCountMismatch { rows: usize, declared: usize },
    /// The water and hexadecane runs describe different numbers of atoms.
    RunLengthMismatch { water: usize, hexadecane: usize },
    /// The hexadecane total area is zero, so the surface-tension coefficient is undefined.
    ZeroArea,
}

impl fmt::Display for SolvError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolvError::AmsolFailed { message } => {
                write!(f, "AMSOL did not complete the calculation: {message}")
            }
            SolvError::NoNameLine => write!(f, "no `<name> <atom count>` line in the AMSOL output"),
            SolvError::NoTotalLine => write!(f, "no `Total:` line in the AMSOL output"),
            SolvError::AtomCountMismatch { rows, declared } => write!(
                f,
                "per-atom table has {rows} rows but the output declares {declared} atoms"
            ),
            SolvError::RunLengthMismatch { water, hexadecane } => write!(
                f,
                "water run has {water} atoms, hexadecane run has {hexadecane}"
            ),
            SolvError::ZeroArea => {
                write!(f, "hexadecane total area is zero")
            }
        }
    }
}

impl std::error::Error for SolvError {}

/// Parse one AMSOL7.1 output file.
///
/// The name line is the first line of exactly two fields whose second field is an integer. Per-atom
/// rows are lines of exactly nine fields whose first field is an integer and whose fifth field is
/// not `*`. The totals line is the first line of exactly six fields beginning `Total:`.
pub fn parse(text: &str) -> Result<AmsolOutput, SolvError> {
    if text.contains(AMSOL_FAILURE_MARKER) {
        // carry AMSOL's own diagnostic rather than reporting a missing section downstream
        let message = text
            .lines()
            .map(str::trim)
            .filter(|l| {
                !l.is_empty() && (l.ends_with('.') || l.starts_with("THE ")) && !l.contains("***")
            })
            .collect::<Vec<_>>()
            .join(" ");
        return Err(SolvError::AmsolFailed { message });
    }
    let mut name: Option<(String, usize)> = None;
    let mut atoms = Vec::new();
    let mut totals: Option<Totals> = None;

    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();

        if name.is_none() && f.len() == 2 {
            if let Ok(n) = f[1].parse::<usize>() {
                name = Some((f[0].to_string(), n));
                continue;
            }
        }

        if f.len() == 9 && f[4] != "*" {
            if let Ok(index) = f[0].parse::<usize>() {
                atoms.push(AtomRow {
                    index,
                    symbol: f[1].to_string(),
                    charge: num(f[2]),
                    g_p: num(f[3]),
                    area: num(f[4]),
                    sigma: num(f[5]),
                    g_cds: num(f[6]),
                    subtotal: num(f[7]),
                    m_value: num(f[8]),
                });
                continue;
            }
        }

        if totals.is_none() && f.len() == 6 && f[0] == "Total:" {
            totals = Some(Totals {
                charge: num(f[1]),
                g_p: num(f[2]),
                area: num(f[3]),
                g_cds: num(f[4]),
                subtotal: num(f[5]),
            });
        }
    }

    let (name, n_atoms) = name.ok_or(SolvError::NoNameLine)?;
    let totals = totals.ok_or(SolvError::NoTotalLine)?;
    if atoms.len() != n_atoms {
        return Err(SolvError::AtomCountMismatch {
            rows: atoms.len(),
            declared: n_atoms,
        });
    }
    Ok(AmsolOutput {
        name,
        n_atoms,
        atoms,
        totals,
    })
}

#[inline]
fn num(s: &str) -> f64 {
    s.parse().unwrap_or(0.0)
}

/// One atom's row in the `.solv` file.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolvAtom {
    /// CM2 charge from the hexadecane run.
    pub charge: f64,
    /// Polar contribution, water minus hexadecane.
    pub diff_pol: f64,
    /// Area from the hexadecane run.
    pub area: f64,
    /// Apolar contribution, water minus hexadecane, less the surface-tension term.
    pub diff_apol: f64,
    /// `diff_pol + diff_apol`.
    pub diff_atomic_solv: f64,
}

/// The contents of a `.solv` file.
#[derive(Clone, Debug, PartialEq)]
pub struct SolvFile {
    pub name: String,
    /// Net charge, from the water run's totals.
    pub charge: f64,
    /// Total area, from the water run's totals.
    pub area: f64,
    pub tot_diff_pol: f64,
    pub tot_diff_apol: f64,
    pub tot_diff_pol_plus_apol: f64,
    pub atoms: Vec<SolvAtom>,
}

/// Form the water-minus-hexadecane differences for one molecule.
///
/// The apolar difference subtracts a surface-tension term present only in the hexadecane run:
/// `cs_coeff = (total G_CDS − Σ atomic G_CDS) / total area`, applied per atom as `cs_coeff × area`.
pub fn combine(water: &AmsolOutput, hexadecane: &AmsolOutput) -> Result<SolvFile, SolvError> {
    if water.atoms.len() != hexadecane.atoms.len() {
        return Err(SolvError::RunLengthMismatch {
            water: water.atoms.len(),
            hexadecane: hexadecane.atoms.len(),
        });
    }
    if hexadecane.totals.area == 0.0 {
        return Err(SolvError::ZeroArea);
    }

    let sum_apol_hex: f64 = hexadecane.atoms.iter().map(|a| a.g_cds).sum();
    let cs_coeff = (hexadecane.totals.g_cds - sum_apol_hex) / hexadecane.totals.area;

    let mut atoms = Vec::with_capacity(water.atoms.len());
    let (mut tot_diff_pol, mut tot_diff_apol, mut tot_diff_pol_plus_apol) = (0.0, 0.0, 0.0);

    for (w, h) in water.atoms.iter().zip(&hexadecane.atoms) {
        let diff_pol = w.g_p - h.g_p;
        let diff_apol = w.g_cds - (h.g_cds + cs_coeff * h.area);
        let diff_atomic_solv = diff_pol + diff_apol;

        tot_diff_pol += diff_pol;
        tot_diff_apol += diff_apol;
        tot_diff_pol_plus_apol += diff_atomic_solv;

        atoms.push(SolvAtom {
            charge: h.charge,
            diff_pol,
            area: h.area,
            diff_apol,
            diff_atomic_solv,
        });
    }

    Ok(SolvFile {
        name: water.name.clone(),
        charge: water.totals.charge,
        area: water.totals.area,
        tot_diff_pol,
        tot_diff_apol,
        tot_diff_pol_plus_apol,
        atoms,
    })
}

impl SolvFile {
    /// Render in the `.solv` layout: a header line of
    /// `%s %3d %4.1f %8.2f %8.2f %8.2f %8.2f`, then one `%8.4f%8.2f%7.2f%8.2f%8.2f` row per atom.
    pub fn render(&self) -> String {
        let mut s = format!(
            "{} {:3} {:4.1} {:8.2} {:8.2} {:8.2} {:8.2}\n",
            self.name,
            self.atoms.len(),
            self.charge,
            self.tot_diff_pol,
            self.area,
            self.tot_diff_apol,
            self.tot_diff_pol_plus_apol,
        );
        for a in &self.atoms {
            s.push_str(&format!(
                "{:8.4}{:8.2}{:7.2}{:8.2}{:8.2}\n",
                a.charge, a.diff_pol, a.area, a.diff_apol, a.diff_atomic_solv
            ));
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal AMSOL output carrying the three sections [`parse`] reads. Synthetic: it exercises
    /// the parser's structure, and is deliberately NOT a reference for any computed value.
    const SYNTHETIC: &str = "\
 mol0001 3

         Atom NO.   Type          Charge        No. of electrons
    1    C   -0.200   0.100  20.00   0.010   0.300   0.400   1.0
    2    H    0.100  -0.050  10.00   0.020  -0.100  -0.150   1.0
    3    H    0.100  -0.050  10.00   0.020  -0.100  -0.150   1.0
 Total:   0.00   0.00  40.00   0.10   0.10
";

    fn synthetic() -> AmsolOutput {
        parse(SYNTHETIC).expect("parsed")
    }

    #[test]
    fn parses_name_atom_table_and_totals() {
        let o = synthetic();
        assert_eq!(o.name, "mol0001");
        assert_eq!(o.n_atoms, 3);
        assert_eq!(o.atoms.len(), 3);
        assert_eq!(o.atoms[0].symbol, "C");
        assert_eq!(o.atoms[0].charge, -0.2);
        assert_eq!(o.totals.area, 40.0);
    }

    #[test]
    fn per_atom_charge_and_area_come_from_the_hexadecane_run() {
        let w = synthetic();
        let mut h = synthetic();
        for a in &mut h.atoms {
            a.charge += 1.0;
            a.area += 5.0;
        }
        let s = combine(&w, &h).expect("combined");
        for (row, hx) in s.atoms.iter().zip(&h.atoms) {
            assert_eq!(row.charge, hx.charge);
            assert_eq!(row.area, hx.area);
        }
        // the header's charge and area come from the water totals
        assert_eq!(s.charge, w.totals.charge);
        assert_eq!(s.area, w.totals.area);
    }

    #[test]
    fn totals_are_the_sums_of_the_per_atom_columns() {
        let (w, h) = (synthetic(), synthetic());
        let s = combine(&w, &h).expect("combined");
        let sum_pol: f64 = s.atoms.iter().map(|a| a.diff_pol).sum();
        let sum_apol: f64 = s.atoms.iter().map(|a| a.diff_apol).sum();
        assert!((s.tot_diff_pol - sum_pol).abs() < 1e-12);
        assert!((s.tot_diff_apol - sum_apol).abs() < 1e-12);
        assert!((s.tot_diff_pol_plus_apol - (sum_pol + sum_apol)).abs() < 1e-12);
    }

    #[test]
    fn render_uses_the_solv_column_layout() {
        let (w, h) = (synthetic(), synthetic());
        let s = combine(&w, &h).expect("combined");
        let rendered = s.render();
        let lines: Vec<&str> = rendered.lines().collect();
        assert_eq!(lines.len(), 1 + s.atoms.len());
        assert!(
            lines[0].starts_with("mol0001   3 "),
            "header: {:?}",
            lines[0]
        );
        // each atom row is the fixed 8+8+7+8+8 layout
        assert!(
            lines[1..].iter().all(|l| l.len() == 39),
            "row width: {:?}",
            lines[1]
        );
    }

    #[test]
    fn mismatched_runs_are_rejected() {
        let w = synthetic();
        let mut h = synthetic();
        h.atoms.pop();
        assert_eq!(
            combine(&w, &h),
            Err(SolvError::RunLengthMismatch {
                water: 3,
                hexadecane: 2
            })
        );
    }

    #[test]
    fn missing_sections_are_reported() {
        assert_eq!(parse(""), Err(SolvError::NoNameLine));
        assert_eq!(parse("mol 5\n"), Err(SolvError::NoTotalLine));
    }

    #[test]
    fn declared_count_must_match_the_table() {
        let doctored = SYNTHETIC.replace("mol0001 3", "mol0001 4");
        assert_eq!(
            parse(&doctored),
            Err(SolvError::AtomCountMismatch {
                rows: 3,
                declared: 4
            })
        );
    }

    #[test]
    fn amsol_size_limit_is_reported_as_such() {
        // AMSOL exits 0 when it hits a compiled-in limit; only the output text says otherwise.
        let text = "\n ***WARNING***\nTHE MAXIMUM NUMBER OF HYDROGEN ATOMS HAS BEEN EXCEEDED.\n\
                    The submitted job was not completed successfully.\n";
        match parse(text) {
            Err(SolvError::AmsolFailed { message }) => {
                assert!(
                    message.contains("MAXIMUM NUMBER OF HYDROGEN ATOMS"),
                    "{message}"
                )
            }
            other => panic!("expected AmsolFailed, got {other:?}"),
        }
    }
}
