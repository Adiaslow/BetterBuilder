//! Assigning experimental torsion preferences.
//!
//! Reproduces `getExperimentalTorsions`. Patterns are tried in library order and the first one to
//! match a given central bond claims it — that ordering is how "most specific wins" is encoded, so
//! the library must be read in file order and matched in that order.
//!
//! Bonds of small rings fused to another ring by more than one bond are excluded: a torsion there
//! is determined by the ring closure, not by the pattern.

use crate::smarts_ast::Pattern;
use crate::smarts_match::{find_matches_with, pattern_graph, Matcher, Perceived};
use crate::torsion_lib::TorsionPattern;
use crate::vf2;

/// Rings smaller than this are not macrocycles, for the fused-ring exclusion.
const MIN_MACROCYCLE_SIZE: usize = 9;

/// One assigned torsion.
#[derive(Debug, Clone, PartialEq)]
pub struct Torsion {
    pub atoms: [usize; 4],
    pub signs: [i8; 6],
    pub v: [f64; 6],
}

/// Ring membership expressed over bonds: one entry per ring, listing its bond indices.
fn bond_rings(mol: &Perceived) -> Vec<Vec<usize>> {
    let mut index = std::collections::HashMap::new();
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        index.insert((a.min(b), a.max(b)), bi);
    }
    mol.rings
        .iter()
        .map(|ring| {
            let set: std::collections::HashSet<usize> = ring.iter().copied().collect();
            mol.bonds
                .iter()
                .enumerate()
                .filter(|(_, &(a, b))| set.contains(&a) && set.contains(&b))
                .map(|(bi, _)| bi)
                .filter(|bi| mol.bond_in_ring.get(*bi).copied().unwrap_or(false))
                .collect()
        })
        .collect()
}

/// Bonds whose torsion is fixed by a fused small ring rather than by a pattern.
fn excluded_bonds(brings: &[Vec<usize>], n_bonds: usize) -> Vec<bool> {
    let mut excluded = vec![false; n_bonds];
    for (i, ri) in brings.iter().enumerate() {
        let si: std::collections::HashSet<usize> = ri.iter().copied().collect();
        for rj in brings.iter().skip(i + 1) {
            if ri.len() >= MIN_MACROCYCLE_SIZE && rj.len() >= MIN_MACROCYCLE_SIZE {
                continue;
            }
            let shared = rj.iter().filter(|b| si.contains(b)).count();
            if shared > 1 {
                if ri.len() < MIN_MACROCYCLE_SIZE {
                    for &b in ri {
                        excluded[b] = true;
                    }
                }
                if rj.len() < MIN_MACROCYCLE_SIZE {
                    for &b in rj {
                        excluded[b] = true;
                    }
                }
            }
        }
    }
    excluded
}

/// A library pattern parsed and with its query graph built, ready to match against any molecule.
pub struct Compiled {
    pub pattern: Pattern,
    pub graph: vf2::Graph,
    /// the torsion quadruple by map number, or `None` if the pattern carries no `:1..:4`
    pub quad: Option<[usize; 4]>,
}

/// Assign torsions, given the library and its compiled patterns in the same order.
pub fn assign(mol: &Perceived, lib: &[TorsionPattern], asts: &[Compiled]) -> (Vec<Torsion>, Vec<bool>) {
    let matcher = Matcher::new(mol.view());
    let n_bonds = mol.bonds.len();
    let brings = bond_rings(mol);
    let excluded = excluded_bonds(&brings, n_bonds);
    let mut rings_per_bond = vec![0usize; n_bonds];
    for r in &brings {
        for &b in r {
            rings_per_bond[b] += 1;
        }
    }
    let mut bond_of: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        bond_of.insert((a.min(b), a.max(b)), bi);
    }

    let mut done = vec![false; n_bonds];
    let mut out = Vec::new();

    for (param, comp) in lib.iter().zip(asts) {
        let Some(quad) = comp.quad else { continue };
        for m in find_matches_with(&comp.pattern, &comp.graph, &matcher) {
            let a: Vec<usize> = quad.iter().map(|&q| m[q]).collect();
            let key = (a[1].min(a[2]), a[1].max(a[2]));
            let Some(&bid) = bond_of.get(&key) else { continue };
            if excluded[bid] || rings_per_bond[bid] > 3 {
                done[bid] = true;
            }
            if done[bid] {
                continue;
            }
            out.push(Torsion {
                atoms: [a[0], a[1], a[2], a[3]],
                signs: param.signs,
                v: param.v,
            });
            done[bid] = true;
        }
    }
    (out, done)
}

/// Parse the library once and build each pattern's query graph; reused across every molecule.
pub fn compile(lib: &[TorsionPattern]) -> Vec<Compiled> {
    lib.iter()
        .map(|p| {
            let pattern = crate::smarts::parse(&p.smarts).expect("library pattern parses");
            let graph = pattern_graph(&pattern);
            let quad = pattern.torsion_quad();
            Compiled { pattern, graph, quad }
        })
        .collect()
}

/// One UFF inversion centre, as `getExperimentalTorsions` records it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Improper {
    /// `[neighbour0, centre, neighbour1, neighbour2]`
    pub atoms: [usize; 4],
    pub atomic_num: u8,
    pub is_c_bound_to_sp2_o: bool,
}

/// Put a ring's atoms in cyclic order. The SSSR stores them sorted, but the planarity terms below
/// walk consecutive atoms, so the cycle has to be recovered from adjacency.
fn cyclic_order(ring: &[usize], adj: &[Vec<usize>]) -> Option<Vec<usize>> {
    let set: std::collections::HashSet<usize> = ring.iter().copied().collect();
    let mut path = vec![*ring.first()?];
    let mut seen: std::collections::HashSet<usize> = path.iter().copied().collect();
    while path.len() < ring.len() {
        let last = *path.last()?;
        let next = adj[last]
            .iter()
            .copied()
            .find(|n| set.contains(n) && !seen.contains(n))?;
        seen.insert(next);
        path.push(next);
    }
    // the walk must close back on itself
    if adj[*path.last()?].contains(&path[0]) {
        Some(path)
    } else {
        None
    }
}

/// Basic-knowledge terms: planarity around SP2 ring bonds, and the inversion centres.
///
/// Runs after the library patterns, sharing their `done` set — a bond already claimed by a pattern
/// is not given a planarity term.
pub fn basic_knowledge(
    mol: &Perceived,
    done: &mut [bool],
) -> (Vec<Torsion>, Vec<Improper>) {
    use crate::hybrid::Hybridization::Sp2;
    let n = mol.atomic_numbers.len();
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in &mol.bonds {
        adj[a].push(b);
        adj[b].push(a);
    }
    let mut bond_of: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        bond_of.insert((a.min(b), a.max(b)), bi);
    }

    // inversion centres: SP2 C/N/O with three neighbours
    let mut impropers = Vec::new();
    for c in 0..n {
        let z = mol.atomic_numbers[c];
        if !matches!(z, 6 | 7 | 8) || mol.hybridization[c] != Sp2 || adj[c].len() != 3 {
            continue;
        }
        let nb = &adj[c];
        let is_c_bound_to_sp2_o = z == 6
            && nb
                .iter()
                .any(|&x| mol.atomic_numbers[x] == 8 && mol.hybridization[x] == Sp2);
        impropers.push(Improper {
            atoms: [nb[0], c, nb[1], nb[2]],
            atomic_num: z,
            is_c_bound_to_sp2_o,
        });
    }

    // planarity across SP2 bonds in rings of size 4..=6
    let mut torsions = Vec::new();
    for ring in &mol.rings {
        let size = ring.len();
        if !(4..=6).contains(&size) {
            continue;
        }
        let Some(cyc) = cyclic_order(ring, &adj) else { continue };
        for i in 0..size {
            let a = [
                cyc[i],
                cyc[(i + 1) % size],
                cyc[(i + 2) % size],
                cyc[(i + 3) % size],
            ];
            let key = (a[1].min(a[2]), a[1].max(a[2]));
            let Some(&bid) = bond_of.get(&key) else { continue };
            if done[bid] || a.iter().any(|&x| mol.hybridization[x] != Sp2) {
                continue;
            }
            done[bid] = true;
            torsions.push(Torsion {
                atoms: a,
                signs: [1, -1, 1, 1, 1, 1],
                v: [0.0, 100.0, 0.0, 0.0, 0.0, 0.0],
            });
        }
    }

    (torsions, impropers)
}
