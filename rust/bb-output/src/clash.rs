//! Clash detection — a faithful port of `clash.py`'s `decideDistanceRules` + the per-set loop in
//! `_identifyClashSetnums`. A set (input conformation) is *broken* when any atom pair matching a rule
//! violates its distance constraint. The pipeline's clashfile (vendored below) holds a single rule,
//! `min 2 1 H H 1.8`: for hydrogen pairs more than 2 bonds apart, a distance under 1.8 Å is a clash.
//! Bond separation is the shortest bond path (BFS over the bond graph), matching `bondsBetweenActual`.

use crate::mol::TypedMol;
use std::collections::VecDeque;

/// The pipeline's clash rules, verbatim (`$DOCKBASE/ligand/mol2db2/clashfile.txt`).
const CLASHFILE: &str = include_str!("clashfile.txt");

struct Rule {
    is_min: bool,   // "min" vs "max"
    bondc: i64,     // bond-count threshold
    cmp: i32,       // -1 lt, 0 eq, 1 gt (how bondsBetween must relate to bondc)
    type_a: String, // atom-type prefix or "*"
    type_b: String,
    dist2: f64, // squared distance constraint
}

fn parse_rules(text: &str) -> Vec<Rule> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            let d: f64 = t[5].parse().unwrap();
            Rule {
                is_min: t[0] == "min",
                bondc: t[1].parse().unwrap(),
                cmp: t[2].parse().unwrap(),
                type_a: t[3].to_string(),
                type_b: t[4].to_string(),
                dist2: d * d,
            }
        })
        .collect()
}

fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    let (dx, dy, dz) = (b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    dx * dx + dy * dy + dz * dz
}

/// Shortest bond-path distances from `start` to every atom (BFS; `-1` if unreachable), reproducing
/// `bondsBetweenActual` for an unweighted bond graph.
fn bond_dist_from(mol2: &TypedMol, start: usize) -> Vec<i64> {
    let n = mol2.atom_bonds.len();
    let mut d = vec![-1i64; n];
    d[start] = 0;
    let mut q = VecDeque::from([start]);
    while let Some(a) = q.pop_front() {
        for (b, _) in &mol2.atom_bonds[a] {
            if d[*b] == -1 {
                d[*b] = d[a] + 1;
                q.push_back(*b);
            }
        }
    }
    d
}

fn matches(cmp: i32, bonds: i64, bondc: i64) -> bool {
    match cmp {
        -1 => bonds < bondc,
        1 => bonds > bondc,
        _ => bonds == bondc,
    }
}

/// One atom pair a rule constrains, with the distance test it demands. Everything but the coordinate
/// comparison is geometry-independent, so these are built once ([`clash_pairs`]) and reused for every
/// conformer.
struct ClashPair {
    a: usize,
    b: usize,
    is_min: bool,
    dist2: f64,
}

/// All atom pairs the rules constrain — the geometry-independent part of `decideDistanceRules` (atom
/// types + bond distances) evaluated once. Mirrors the pair enumeration the old per-conformer
/// `conf_clashes` did inline: same rule order, same three loops, same `matches` filter. Only the
/// `typeA == typeB` branch is exercised in practice (the sole rule is H–H); the others are ported for
/// fidelity.
fn clash_pairs(mol2: &TypedMol, rules: &[Rule], bond_dist: &[Vec<i64>]) -> Vec<ClashPair> {
    let natoms = mol2.atom_type.len();
    let of_type = |ty: &str| -> Vec<usize> {
        (0..natoms).filter(|&a| ty == "*" || mol2.atom_type[a].starts_with(ty)).collect()
    };
    let mut pairs = Vec::new();
    let consider = |a: usize, b: usize, rule: &Rule, pairs: &mut Vec<ClashPair>| {
        if matches(rule.cmp, bond_dist[a][b], rule.bondc) {
            pairs.push(ClashPair { a, b, is_min: rule.is_min, dist2: rule.dist2 });
        }
    };
    for rule in rules {
        if rule.type_a == rule.type_b {
            let atoms = of_type(&rule.type_a);
            for i in 0..atoms.len() {
                for j in (i + 1)..atoms.len() {
                    consider(atoms[i], atoms[j], rule, &mut pairs);
                }
            }
        } else {
            // typeA != typeB: check the defined type against itself then against all others
            let (a_star, b_star) = (rule.type_a == "*", rule.type_b == "*");
            let (group_a, group_b): (Vec<usize>, Vec<usize>) = if a_star ^ b_star {
                let def = if a_star { &rule.type_b } else { &rule.type_a };
                let defs = of_type(def);
                let others: Vec<usize> = (0..natoms).filter(|&a| !mol2.atom_type[a].starts_with(def)).collect();
                (defs, others)
            } else {
                (of_type(&rule.type_a), of_type(&rule.type_a))
            };
            // defined-vs-defined
            for i in 0..group_a.len() {
                for j in (i + 1)..group_a.len() {
                    consider(group_a[i], group_a[j], rule, &mut pairs);
                }
            }
            // defined-vs-others
            for &a in &group_a {
                for &b in &group_b {
                    consider(a, b, rule, &mut pairs);
                }
            }
        }
    }
    pairs
}

/// True iff conformation `xyz` violates any precomputed clash pair.
fn conf_clashes(xyz: &[[f64; 3]], pairs: &[ClashPair]) -> bool {
    pairs.iter().any(|p| {
        let dd = dist2(xyz[p.a], xyz[p.b]);
        (p.is_min && dd < p.dist2) || (!p.is_min && dd > p.dist2)
    })
}

/// The set indices (input conformations) that are broken, in ascending order.
pub fn broken_sets(mol2: &TypedMol) -> Vec<usize> {
    let rules = parse_rules(CLASHFILE);
    let bond_dist: Vec<Vec<i64>> = (0..mol2.atom_bonds.len()).map(|a| bond_dist_from(mol2, a)).collect();
    let pairs = clash_pairs(mol2, &rules, &bond_dist); // geometry-independent, built once
    (0..mol2.atom_xyz.len())
        .filter(|&s| conf_clashes(&mol2.atom_xyz[s], &pairs))
        .collect()
}
