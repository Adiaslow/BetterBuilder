//! Torsion-strain computation — the runtime half of the port of `TL_Functions.TL_Lookup`.
//!
//! Matching is geometry-independent, so it runs once per molecule ([`match_library`]); the per-
//! conformer energy ([`strain`]) then takes each matched torsion's dihedral and reads its penalty.
//! Two reductions keep one torsion per rotatable bond: first the most-specific rule for identical
//! 4-atom torsions, then, per central bond, specific-over-general and higher-energy. Total strain is
//! the sum of the non-negative per-torsion energies (or `−(#flagged)` if any approximate angle is
//! unobserved); max strain is the largest single-torsion energy.

use crate::library::{library, Method, Rule};
use bb_perceive::smarts;
use bb_perceive::smarts_match::{find_matches_with, pattern_graph, Matcher, Perceived};

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = dot(v, v).sqrt();
    if n == 0.0 {
        v
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Dihedral angle in degrees, (−180, 180], matching `TL_Functions.dihedral` (incl. its negation).
fn dihedral(a1: [f64; 3], a2: [f64; 3], a3: [f64; 3], a4: [f64; 3]) -> f64 {
    let (b1, b2, b3) = (sub(a2, a1), sub(a3, a2), sub(a4, a3));
    let n1 = unit(cross(b1, b2));
    let n2 = unit(cross(b2, b3));
    let m = unit(cross(n1, b2));
    let x = dot(n1, n2);
    let y = dot(m, n2);
    -y.atan2(x) * 180.0 / std::f64::consts::PI
}

/// Signed angular difference `theta_1 − theta_2`, wrapped to (−180, 180], matching `ang_diff`.
fn ang_diff(t1: f64, t2: f64) -> f64 {
    let mut d = (t1 - t2).rem_euclid(360.0);
    if d > 180.0 {
        d -= 360.0;
    }
    d
}

/// One matched torsion: the four atoms (:1..:4 order) and the rule it matched.
pub struct Hit {
    pub atoms: [usize; 4],
    pub rule: &'static Rule,
}

/// All library-rule matches on a molecule (geometry-independent). Reproduces `get_matches_idx`: the
/// first four SMARTS atoms are the torsion; matches spanning more than four atoms, or whose end atoms
/// are hydrogens, are skipped. Rules whose SMARTS bb-perceive cannot parse are skipped (surfaced by the
/// gate if any matter).
///
/// All patterns are compiled up front so every one (and its recursive `$(...)` sub-patterns) stays
/// alive with a stable, distinct address for the whole match: `Matcher`'s recursive-match cache is
/// keyed by sub-pattern pointer, so parsing patterns one-at-a-time (freeing between) would let a reused
/// address return a previous rule's cached result.
pub fn match_library(p: &Perceived) -> Vec<Hit> {
    let compiled: Vec<(bb_perceive::smarts_ast::Pattern, bb_perceive::vf2::Graph, &'static Rule)> = library()
        .iter()
        .filter_map(|rule| {
            let pat = smarts::parse(&rule.smarts).ok()?;
            let g = pattern_graph(&pat);
            Some((pat, g, rule))
        })
        .collect();
    let matcher = Matcher::new(p.view());
    let mut hits = Vec::new();
    for (pat, g, rule) in &compiled {
        for m in find_matches_with(pat, g, &matcher) {
            if m.len() > 4 {
                continue;
            }
            if p.atomic_numbers[m[0]] == 1 || p.atomic_numbers[m[3]] == 1 {
                continue;
            }
            hits.push(Hit { atoms: [m[0], m[1], m[2], m[3]], rule });
        }
    }
    hits
}

fn energy(method: &Method, theta: f64) -> (f64, bool) {
    match method {
        Method::Exact { energy, lower: _, upper: _ } => {
            let bin_num = (theta / 10.0).ceil() as i64 + 17;
            let bn = bin_num as usize;
            let prev = ((bin_num + 35) % 36) as usize;
            let e = (energy[bn] - energy[prev]) / 10.0 * (theta - (bin_num - 17) as f64 * 10.0) + energy[bn];
            (e, false)
        }
        Method::Approximate { angles } => {
            for a in angles {
                let delta = ang_diff(theta, a.theta_0);
                if delta.abs() <= a.tolerance2 {
                    return (a.beta_1 * delta.powi(2) + a.beta_2 * delta.powi(4), false);
                }
            }
            (100.0, true) // not observed
        }
    }
}

/// A reduced per-bond torsion: canonicalized atoms, energy, source rule.
pub struct BondInfo {
    pub atoms: [usize; 4], // direction-canonicalized (atoms[1] < atoms[2])
    pub energy: f64,
    pub is_general: bool,
    pub index: usize,
    pub not_observed: bool,
    pub smarts: &'static str,
}

/// The two-stage reduction over a conformer's matched torsions — one entry per rotatable bond.
pub fn reduce(hits: &[Hit], coords: &[[f64; 3]]) -> Vec<BondInfo> {
    let mut bond_info: Vec<BondInfo> = Vec::new();
    for h in hits {
        let theta = dihedral(coords[h.atoms[0]], coords[h.atoms[1]], coords[h.atoms[2]], coords[h.atoms[3]]);
        let (e, not_observed) = energy(&h.rule.method, theta);
        let mut atoms = h.atoms;
        if atoms[1] > atoms[2] {
            atoms.reverse(); // canonicalize bond direction, as the source does
        }
        bond_info.push(BondInfo {
            atoms,
            energy: e,
            is_general: h.rule.is_general,
            index: h.rule.index,
            not_observed,
            smarts: &h.rule.smarts,
        });
    }

    // reduction 1: one entry per identical 4-atom torsion, keeping the lowest (most specific) index
    let mut red1: Vec<BondInfo> = Vec::new();
    for b in bond_info {
        match red1.iter_mut().find(|r| r.atoms == b.atoms) {
            Some(r) => {
                if b.index < r.index {
                    *r = b;
                }
            }
            None => red1.push(b),
        }
    }

    // reduction 2: one entry per central bond, specific over general, then higher energy
    let mut red2: Vec<BondInfo> = Vec::new();
    for b in red1 {
        match red2.iter_mut().find(|r| r.atoms[1] == b.atoms[1] && r.atoms[2] == b.atoms[2]) {
            Some(r) => {
                let more_specific = !b.is_general && r.is_general;
                let higher_same = b.is_general == r.is_general && b.energy > r.energy;
                if more_specific || higher_same {
                    *r = b;
                }
            }
            None => red2.push(b),
        }
    }
    red2
}

/// Per-conformer (total, max) strain for a molecule's precomputed matches at the given coordinates.
pub fn strain(hits: &[Hit], coords: &[[f64; 3]]) -> (f64, f64) {
    let red2 = reduce(hits, coords);
    if red2.is_empty() {
        return (0.0, 0.0);
    }
    let flagged = red2.iter().filter(|b| b.not_observed).count();
    let total = if flagged > 0 {
        -(flagged as f64)
    } else {
        red2.iter().map(|b| b.energy).filter(|&e| e >= 0.0).sum()
    };
    let max = red2.iter().map(|b| b.energy).fold(f64::NEG_INFINITY, f64::max);
    (total, max)
}
