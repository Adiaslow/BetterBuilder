//! Cartesian → MOPAC internal coordinates (Z-matrix), the geometry format carried by the AMSOL7.1
//! input files.
//!
//! Each atom past the first is expressed as a distance to one earlier atom, an angle to a second,
//! and a dihedral to a third. Atom order is preserved: line `i` describes input atom `i`, because
//! AMSOL reports its per-atom results in input order and those are matched positionally against a
//! mol2 written in the same order.

use std::fmt;

/// A reference is rejected when the sine of the angle it forms falls below this, leaving the angle
/// or dihedral it defines numerically ill-conditioned.
const MIN_SIN: f64 = 0.05;

#[derive(Debug, PartialEq, Eq)]
pub enum ZmatrixError {
    /// `atomic_numbers` and `coords` imply different atom counts.
    LengthMismatch { atoms: usize, coords: usize },
    /// Fewer than one atom.
    Empty,
    /// An atomic number with no corresponding element.
    UnknownElement { atom: usize, z: u8 },
    /// Two atoms occupy the same position, leaving a bond length of zero.
    CoincidentAtoms { a: usize, b: usize },
}

impl fmt::Display for ZmatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ZmatrixError::LengthMismatch { atoms, coords } => write!(
                f,
                "{atoms} atomic numbers but {coords} coordinates (expected {})",
                atoms * 3
            ),
            ZmatrixError::Empty => write!(f, "no atoms"),
            ZmatrixError::UnknownElement { atom, z } => {
                write!(f, "atom {atom}: {z} is not an atomic number")
            }
            ZmatrixError::CoincidentAtoms { a, b } => {
                write!(f, "atoms {a} and {b} occupy the same position")
            }
        }
    }
}

impl std::error::Error for ZmatrixError {}

/// One Z-matrix line: the element symbol, the three internal coordinates with their optimization
/// flags, and the 1-based indices of the atoms they are measured against (0 where unused).
#[derive(Clone, Debug, PartialEq)]
pub struct ZLine {
    pub symbol: &'static str,
    pub bond_length: f64,
    pub bond_flag: u8,
    /// Degrees.
    pub angle: f64,
    pub angle_flag: u8,
    /// Degrees.
    pub dihedral: f64,
    pub dihedral_flag: u8,
    pub na: usize,
    pub nb: usize,
    pub nc: usize,
}

impl fmt::Display for ZLine {
    /// The column layout AMSOL's input files carry:
    /// `%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:<2} {:10.6} {:2} {:11.6} {:2} {:11.6} {:2} {:5} {:3} {:3}",
            self.symbol,
            self.bond_length,
            self.bond_flag,
            self.angle,
            self.angle_flag,
            self.dihedral,
            self.dihedral_flag,
            self.na,
            self.nb,
            self.nc,
        )
    }
}

#[inline]
fn pt(coords: &[f64], i: usize) -> [f64; 3] {
    [coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]]
}
#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Distance between atoms `i` and `j`.
fn distance(coords: &[f64], i: usize, j: usize) -> f64 {
    norm(sub(pt(coords, i), pt(coords, j)))
}

/// Angle i-j-k in degrees, `j` central.
fn angle_deg(coords: &[f64], i: usize, j: usize, k: usize) -> f64 {
    let (u, v) = (
        sub(pt(coords, i), pt(coords, j)),
        sub(pt(coords, k), pt(coords, j)),
    );
    let m = norm(u) * norm(v);
    if m < 1e-12 {
        return 0.0;
    }
    (dot(u, v) / m).clamp(-1.0, 1.0).acos().to_degrees()
}

/// Sine of the angle i-j-k, used to score how well-conditioned a reference choice is.
fn angle_sin(coords: &[f64], i: usize, j: usize, k: usize) -> f64 {
    angle_deg(coords, i, j, k).to_radians().sin().abs()
}

/// Signed dihedral i-j-k-l in degrees, in `(-180, 180]`, using the same convention as the
/// engine's torsion term.
fn dihedral_deg(coords: &[f64], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let (pi, pj, pk, pl) = (pt(coords, i), pt(coords, j), pt(coords, k), pt(coords, l));
    let (b1, b2, b3) = (sub(pj, pi), sub(pk, pj), sub(pl, pk));
    let (n1, n2) = (cross(b1, b2), cross(b2, b3));
    let b2n = norm(b2);
    if norm(n1) < 1e-12 || norm(n2) < 1e-12 || b2n < 1e-12 {
        return 0.0;
    }
    (dot(cross(n1, n2), b2) / b2n)
        .atan2(dot(n1, n2))
        .to_degrees()
}

/// Candidate reference atoms for atom `i`, ordered: atoms bonded to `anchor` first (ascending
/// index), then every other already-placed atom by increasing distance from `i`.
fn candidates(i: usize, anchor: usize, coords: &[f64], adj: &[Vec<usize>]) -> Vec<usize> {
    let mut bonded: Vec<usize> = adj[anchor].iter().copied().filter(|&c| c < i).collect();
    bonded.sort_unstable();
    let mut rest: Vec<usize> = (0..i).filter(|c| !bonded.contains(c)).collect();
    rest.sort_by(|&a, &b| {
        distance(coords, i, a)
            .partial_cmp(&distance(coords, i, b))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    bonded.extend(rest);
    bonded
}

/// The first candidate scoring at least [`MIN_SIN`], or the highest-scoring one when none does.
fn pick(cands: &[usize], exclude: &[usize], score: impl Fn(usize) -> f64) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for &c in cands.iter().filter(|c| !exclude.contains(c)) {
        let s = score(c);
        if s >= MIN_SIN {
            return Some(c);
        }
        if best.is_none_or(|(_, bs)| s > bs) {
            best = Some((c, s));
        }
    }
    best.map(|(c, _)| c)
}

/// Convert Cartesian coordinates to Z-matrix lines, one per atom, in input order.
///
/// `atomic_numbers` is one entry per atom, `coords` is `3 × n_atoms` row-major xyz in Å, and
/// `bonds` lists atom index pairs. References are drawn from the bond graph where it supplies a
/// well-conditioned choice and from the nearest already-placed atoms otherwise; `bonds` may be
/// empty, in which case every reference comes from proximity.
///
/// The optimization flags reproduce the AMSOL input files' fixed pattern: the coordinates that are
/// undefined for the first three atoms are flagged 0, every other coordinate 1.
pub fn to_zmatrix(
    atomic_numbers: &[u8],
    coords: &[f64],
    bonds: &[[u32; 2]],
) -> Result<Vec<ZLine>, ZmatrixError> {
    let n = atomic_numbers.len();
    if n == 0 {
        return Err(ZmatrixError::Empty);
    }
    if coords.len() != n * 3 {
        return Err(ZmatrixError::LengthMismatch {
            atoms: n,
            coords: coords.len(),
        });
    }
    let mut symbols = Vec::with_capacity(n);
    for (atom, &z) in atomic_numbers.iter().enumerate() {
        let element = (z as usize)
            .checked_sub(1)
            .and_then(|i| mendeleev::ALL_ELEMENTS.get(i))
            .ok_or(ZmatrixError::UnknownElement { atom, z })?;
        symbols.push(element.symbol());
    }

    let mut adj = vec![Vec::new(); n];
    for b in bonds {
        let (i, j) = (b[0] as usize, b[1] as usize);
        if i < n && j < n && i != j {
            adj[i].push(j);
            adj[j].push(i);
        }
    }

    let mut lines = Vec::with_capacity(n);
    for i in 0..n {
        if i == 0 {
            lines.push(ZLine {
                symbol: symbols[0],
                bond_length: 0.0,
                bond_flag: 0,
                angle: 0.0,
                angle_flag: 0,
                dihedral: 0.0,
                dihedral_flag: 0,
                na: 0,
                nb: 0,
                nc: 0,
            });
            continue;
        }

        // Distance reference: a bonded already-placed atom, else the nearest already-placed atom.
        let cands_na = candidates(i, i, coords, &adj);
        let na = pick(&cands_na, &[], |_| 1.0).expect("at least one earlier atom exists");
        let r = distance(coords, i, na);
        if r < 1e-9 {
            return Err(ZmatrixError::CoincidentAtoms { a: i, b: na });
        }

        // Angle reference: prefer an atom bonded to `na` that is not collinear with i-na.
        let nb = if i >= 2 {
            pick(&candidates(i, na, coords, &adj), &[na], |c| {
                angle_sin(coords, i, na, c)
            })
        } else {
            None
        };
        let a = nb.map_or(0.0, |nb| angle_deg(coords, i, na, nb));

        // Dihedral reference: prefer an atom bonded to `nb` that is not collinear with na-nb.
        let nc = match (i >= 3, nb) {
            (true, Some(nb)) => pick(&candidates(i, nb, coords, &adj), &[na, nb], |c| {
                angle_sin(coords, na, nb, c)
            }),
            _ => None,
        };
        let d = match (nb, nc) {
            (Some(nb), Some(nc)) => dihedral_deg(coords, i, na, nb, nc),
            _ => 0.0,
        };

        lines.push(ZLine {
            symbol: symbols[i],
            bond_length: r,
            bond_flag: 1,
            angle: a,
            angle_flag: u8::from(nb.is_some()),
            dihedral: d,
            dihedral_flag: u8::from(nc.is_some()),
            na: na + 1,
            nb: nb.map_or(0, |x| x + 1),
            nc: nc.map_or(0, |x| x + 1),
        });
    }
    Ok(lines)
}

/// Render Z-matrix lines as the newline-terminated block the AMSOL input files carry.
pub fn render(lines: &[ZLine]) -> String {
    let mut s = String::new();
    for l in lines {
        s.push_str(&l.to_string());
        s.push('\n');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Rebuild Cartesian coordinates from a Z-matrix (natural extension reference frame), so a
    /// round trip can be checked against the input geometry.
    fn reconstruct(lines: &[ZLine]) -> Vec<f64> {
        let mut p: Vec<[f64; 3]> = Vec::with_capacity(lines.len());
        for (i, l) in lines.iter().enumerate() {
            let q = match i {
                0 => [0.0, 0.0, 0.0],
                1 => [l.bond_length, 0.0, 0.0],
                _ if l.nc == 0 => {
                    // place in the xy-plane using distance + angle only
                    let a = p[l.na - 1];
                    let b = p[l.nb - 1];
                    let ab = sub(b, a);
                    let u = {
                        let n = norm(ab);
                        [ab[0] / n, ab[1] / n, ab[2] / n]
                    };
                    // any vector perpendicular to u
                    let seed = if u[0].abs() < 0.9 {
                        [1.0, 0.0, 0.0]
                    } else {
                        [0.0, 1.0, 0.0]
                    };
                    let mut w = cross(u, seed);
                    let wn = norm(w);
                    w = [w[0] / wn, w[1] / wn, w[2] / wn];
                    let th = l.angle.to_radians();
                    [
                        a[0] + l.bond_length * (u[0] * th.cos() + w[0] * th.sin()),
                        a[1] + l.bond_length * (u[1] * th.cos() + w[1] * th.sin()),
                        a[2] + l.bond_length * (u[2] * th.cos() + w[2] * th.sin()),
                    ]
                }
                _ => {
                    let (a, b, c) = (p[l.na - 1], p[l.nb - 1], p[l.nc - 1]);
                    let (th, ph) = (l.angle.to_radians(), l.dihedral.to_radians());
                    let bc = sub(a, b);
                    let bcn = norm(bc);
                    let u = [bc[0] / bcn, bc[1] / bcn, bc[2] / bcn];
                    let ab = sub(b, c);
                    let mut nrm = cross(ab, u);
                    let nn = norm(nrm);
                    nrm = [nrm[0] / nn, nrm[1] / nn, nrm[2] / nn];
                    let m = cross(nrm, u);
                    let (r, st, ct) = (l.bond_length, th.sin(), th.cos());
                    // local frame: -u along the bond, m/nrm spanning the perpendicular plane
                    [
                        a[0] + r * (-u[0] * ct + m[0] * st * ph.cos() + nrm[0] * st * ph.sin()),
                        a[1] + r * (-u[1] * ct + m[1] * st * ph.cos() + nrm[1] * st * ph.sin()),
                        a[2] + r * (-u[2] * ct + m[2] * st * ph.cos() + nrm[2] * st * ph.sin()),
                    ]
                }
            };
            p.push(q);
        }
        p.into_iter().flatten().collect()
    }

    /// Pairwise distance matrix, invariant to the rigid-body placement a Z-matrix does not fix.
    fn dist_matrix(coords: &[f64]) -> Vec<f64> {
        let n = coords.len() / 3;
        let mut d = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                d.push(distance(coords, i, j));
            }
        }
        d
    }

    /// Staggered n-butane: a known 180° C-C-C-C dihedral.
    fn butane() -> (Vec<u8>, Vec<f64>, Vec<[u32; 2]>) {
        let coords = vec![
            -1.9, 0.0, 0.0, // C0
            -0.6, 0.6, 0.0, // C1
            0.6, -0.6, 0.0, // C2
            1.9, 0.0, 0.0, // C3
        ];
        (vec![6, 6, 6, 6], coords, vec![[0, 1], [1, 2], [2, 3]])
    }

    #[test]
    fn round_trip_preserves_geometry() {
        let (z, coords, bonds) = butane();
        let lines = to_zmatrix(&z, &coords, &bonds).expect("zmatrix built");
        let back = reconstruct(&lines);
        let (a, b) = (dist_matrix(&coords), dist_matrix(&back));
        let max = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert!(max < 1e-6, "distance matrix drifted by {max:.3e}");
    }

    #[test]
    fn round_trip_preserves_geometry_on_a_ring() {
        // Puckered 6-ring plus one substituent, so dihedrals span more than one sign.
        let mut coords = Vec::new();
        let mut z = Vec::new();
        let mut bonds = Vec::new();
        for i in 0..6u32 {
            let th = std::f64::consts::TAU * (i as f64) / 6.0;
            coords.extend_from_slice(&[
                1.5 * th.cos(),
                1.5 * th.sin(),
                if i % 2 == 0 { 0.25 } else { -0.25 },
            ]);
            z.push(6);
            bonds.push([i, (i + 1) % 6]);
        }
        coords.extend_from_slice(&[2.6, 0.9, 0.8]);
        z.push(8);
        bonds.push([0, 6]);

        let lines = to_zmatrix(&z, &coords, &bonds).expect("zmatrix built");
        let back = reconstruct(&lines);
        let (a, b) = (dist_matrix(&coords), dist_matrix(&back));
        let max = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
        assert!(max < 1e-6, "distance matrix drifted by {max:.3e}");
    }

    #[test]
    fn butane_anti_dihedral_is_180_degrees() {
        let (z, coords, bonds) = butane();
        let lines = to_zmatrix(&z, &coords, &bonds).expect("zmatrix built");
        // atom 3 is referenced to 2-1-0, so its dihedral is the C-C-C-C torsion
        assert_eq!((lines[3].na, lines[3].nb, lines[3].nc), (3, 2, 1));
        assert!(
            (lines[3].dihedral.abs() - 180.0).abs() < 1e-6,
            "expected |dihedral| 180, got {}",
            lines[3].dihedral
        );
    }

    #[test]
    fn leading_atoms_carry_the_fixed_flag_pattern() {
        let (z, coords, bonds) = butane();
        let lines = to_zmatrix(&z, &coords, &bonds).expect("zmatrix built");
        let flags: Vec<_> = lines
            .iter()
            .map(|l| (l.bond_flag, l.angle_flag, l.dihedral_flag))
            .collect();
        assert_eq!(flags[0], (0, 0, 0));
        assert_eq!(flags[1], (1, 0, 0));
        assert_eq!(flags[2], (1, 1, 0));
        assert_eq!(flags[3], (1, 1, 1));
        // unused references are index 0
        assert_eq!((lines[0].na, lines[0].nb, lines[0].nc), (0, 0, 0));
        assert_eq!((lines[1].nb, lines[1].nc), (0, 0));
        assert_eq!(lines[2].nc, 0);
    }

    #[test]
    fn line_layout_matches_the_amsol_input_columns() {
        let l = ZLine {
            symbol: "C",
            bond_length: 1.54,
            bond_flag: 1,
            angle: 109.47,
            angle_flag: 1,
            dihedral: -60.0,
            dihedral_flag: 1,
            na: 3,
            nb: 2,
            nc: 1,
        };
        // Expected string produced by the layout the AMSOL input files use, applied to these
        // fields: `'%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d' % ('C', 1.54, 1, 109.47,
        // 1, -60.0, 1, 3, 2, 1)`.
        assert_eq!(
            l.to_string(),
            "C    1.540000  1  109.470000  1  -60.000000  1     3   2   1"
        );
    }

    #[test]
    fn rejects_unknown_elements_and_bad_lengths() {
        assert_eq!(
            to_zmatrix(&[6, 200], &[0.0; 6], &[]),
            Err(ZmatrixError::UnknownElement { atom: 1, z: 200 })
        );
        assert_eq!(
            to_zmatrix(&[6, 6], &[0.0; 3], &[]),
            Err(ZmatrixError::LengthMismatch {
                atoms: 2,
                coords: 3
            })
        );
        assert_eq!(to_zmatrix(&[], &[], &[]), Err(ZmatrixError::Empty));
    }

    #[test]
    fn rejects_coincident_atoms() {
        let r = to_zmatrix(&[6, 6], &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], &[[0, 1]]);
        assert_eq!(r, Err(ZmatrixError::CoincidentAtoms { a: 1, b: 0 }));
    }
}
