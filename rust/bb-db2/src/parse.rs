//! db2 text → [`Db2Entry`] parsing. A faithful reader of the `T/M/A/B/X/R/C/S/D/E` records, matching
//! `pyDOCK/ligand/db2.py` / `dockrs-io`: 1-based ids kept as-is, `M1` name/protname read at their
//! fixed columns (they may carry embedded spaces), everything else whitespace-split. Sets use the
//! header's `#lines` count to bound the member lines that follow.

use crate::model::{Atom, Bond, Conf, Coord, Db2Entry, Rigid, Set};

/// Parse error with the offending line's 1-based number and a short reason.
#[derive(Debug)]
pub struct ParseError {
    pub line: usize,
    pub reason: String,
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "db2 parse error at line {}: {}", self.line, self.reason)
    }
}
impl std::error::Error for ParseError {}

fn parse_at<T: std::str::FromStr>(t: &[&str], i: usize, line: usize, what: &str) -> Result<T, ParseError> {
    t.get(i)
        .ok_or_else(|| ParseError { line, reason: format!("missing field {i} ({what})") })?
        .parse()
        .map_err(|_| ParseError { line, reason: format!("bad {what}: {:?}", t.get(i)) })
}

fn xyz(t: &[&str], i: usize, line: usize) -> Result<[f64; 3], ParseError> {
    Ok([
        parse_at(t, i, line, "x")?,
        parse_at(t, i + 1, line, "y")?,
        parse_at(t, i + 2, line, "z")?,
    ])
}

/// Parse all stacked entries in a db2 file body.
pub fn parse(text: &str) -> Result<Vec<Db2Entry>, ParseError> {
    let mut entries = Vec::new();
    let mut e = Db2Entry::default();
    let mut m_idx = 0usize;
    let mut in_entry = false;
    let mut members_remaining = 0usize;

    for (li0, line) in text.lines().enumerate() {
        let line_no = li0 + 1;
        let tag = line.as_bytes().first().copied();
        let t: Vec<&str> = line.split_whitespace().collect();
        match tag {
            Some(b'T') => {
                e.color_table
                    .push((parse_at(&t, 1, line_no, "color code")?, t.get(2).unwrap_or(&"").to_string()));
            }
            Some(b'M') => {
                in_entry = true;
                m_idx += 1;
                match m_idx {
                    // M1: name at cols 3-18 (0-based 2..18), protname cols 20-28 (0-based 19..28).
                    1 => {
                        e.name = line.get(2..18).unwrap_or("").trim().to_string();
                        e.protname = line.get(19..28).unwrap_or("").trim().to_string();
                    }
                    2 => {
                        e.total_charge = parse_at(&t, 1, line_no, "total charge")?;
                        e.total_polar_solv = parse_at(&t, 2, line_no, "total polar solv")?;
                        e.total_apolar_solv = parse_at(&t, 3, line_no, "total apolar solv")?;
                        e.total_solv = parse_at(&t, 4, line_no, "total solv")?;
                        e.surface_area = parse_at(&t, 5, line_no, "surface area")?;
                    }
                    3 => e.smiles = line.get(2..).unwrap_or("").trim_end().to_string(),
                    4 => e.longname = line.get(2..).unwrap_or("").trim_end().to_string(),
                    5 => e.m5 = parse_at(&t, 1, line_no, "M5 value")?,
                    _ => e.extra_m.push(line.get(2..).unwrap_or("").to_string()),
                }
            }
            Some(b'A') => e.atoms.push(Atom {
                num: parse_at(&t, 1, line_no, "atom num")?,
                name: t.get(2).unwrap_or(&"").to_string(),
                sybyl: t.get(3).unwrap_or(&"").to_string(),
                vdwtype: parse_at(&t, 4, line_no, "vdwtype")?,
                color: parse_at(&t, 5, line_no, "color")?,
                charge: parse_at(&t, 6, line_no, "charge")?,
                polar_solv: parse_at(&t, 7, line_no, "polar solv")?,
                apolar_solv: parse_at(&t, 8, line_no, "apolar solv")?,
                total_solv: parse_at(&t, 9, line_no, "total solv")?,
                surface: parse_at(&t, 10, line_no, "surface")?,
            }),
            Some(b'B') => e.bonds.push(Bond {
                num: parse_at(&t, 1, line_no, "bond num")?,
                a1: parse_at(&t, 2, line_no, "bond a1")?,
                a2: parse_at(&t, 3, line_no, "bond a2")?,
                btype: t.get(4).unwrap_or(&"").to_string(),
            }),
            Some(b'X') => e.coords.push(Coord {
                num: parse_at(&t, 1, line_no, "coord num")?,
                atom: parse_at(&t, 2, line_no, "coord atom")?,
                conf: parse_at(&t, 3, line_no, "coord conf")?,
                xyz: xyz(&t, 4, line_no)?,
            }),
            Some(b'R') => e.rigid.push(Rigid {
                num: parse_at(&t, 1, line_no, "rigid num")?,
                color: parse_at(&t, 2, line_no, "rigid color")?,
                xyz: xyz(&t, 3, line_no)?,
            }),
            Some(b'C') => e.confs.push(Conf {
                num: parse_at(&t, 1, line_no, "conf num")?,
                start: parse_at(&t, 2, line_no, "conf start")?,
                end: parse_at(&t, 3, line_no, "conf end")?,
            }),
            Some(b'S') => {
                if members_remaining == 0 {
                    // header: setnum #lines #confs broken hydrogens totalStrain maxStrain
                    members_remaining = parse_at(&t, 2, line_no, "set #lines")?;
                    e.sets.push(Set {
                        num: parse_at(&t, 1, line_no, "set num")?,
                        broken: parse_at(&t, 4, line_no, "set broken")?,
                        hydrogens: parse_at(&t, 5, line_no, "set hydrogens")?,
                        total_strain: parse_at(&t, 6, line_no, "set total strain")?,
                        max_strain: parse_at(&t, 7, line_no, "set max strain")?,
                        member_lines: Vec::new(),
                    });
                } else {
                    // member: setnum linenum #confs conf1 conf2 …
                    let nconfs: usize = parse_at(&t, 3, line_no, "member #confs")?;
                    let mut confs = Vec::with_capacity(nconfs);
                    for k in 0..nconfs {
                        confs.push(parse_at(&t, 4 + k, line_no, "member conf")?);
                    }
                    let set = e.sets.last_mut().ok_or_else(|| ParseError {
                        line: line_no,
                        reason: "S member before any S header".into(),
                    })?;
                    set.member_lines.push(confs);
                    members_remaining -= 1;
                }
            }
            Some(b'D') => e.cluster_lines.push(line.to_string()),
            Some(b'E') => {
                entries.push(std::mem::take(&mut e));
                m_idx = 0;
                in_entry = false;
                members_remaining = 0;
            }
            _ => {}
        }
    }
    if in_entry {
        entries.push(e);
    }
    Ok(entries)
}
