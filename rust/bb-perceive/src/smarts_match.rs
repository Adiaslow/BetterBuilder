//! Evaluating SMARTS predicates and matching patterns against a molecule.
//!
//! Every predicate is answered from our own perception — rings, hybridization, hydrogen counts —
//! all of which reproduce RDKit exactly. Nothing is delegated to a second interpretation of what
//! `X3` or `H1` mean. The subgraph isomorphism itself is [`crate::vf2`], which reproduces the
//! traversal order RDKit's own VF2 uses, driven by these predicates as node and edge matchers.

use crate::hybrid::Hybridization;
use crate::smarts_ast::{AtomPred, BondPred, Pattern};
use crate::vf2;
use std::collections::HashMap;

/// Everything a SMARTS predicate can ask about a molecule, post-`AddHs`.
pub struct MolView<'a> {
    pub atomic_numbers: &'a [u8],
    pub charges: &'a [i8],
    pub aromatic: &'a [bool],
    pub hybridization: &'a [Hybridization],
    /// neighbours of each atom
    pub adj: Vec<Vec<usize>>,
    /// bonds as `(a, b)`, and their integer order (aromatic recorded separately)
    pub bonds: &'a [(usize, usize)],
    pub bond_order: &'a [u8],
    pub bond_aromatic: &'a [bool],
    pub bond_in_ring: &'a [bool],
    /// ring sizes each atom belongs to
    pub ring_sizes: Vec<Vec<usize>>,
    /// number of rings each atom belongs to
    pub ring_count: Vec<usize>,
    /// number of ring bonds at each atom
    pub ring_bonds_at: Vec<usize>,
}

impl MolView<'_> {
    fn total_h(&self, i: usize) -> usize {
        // post-AddHs every hydrogen is a real neighbour
        self.adj[i]
            .iter()
            .filter(|&&n| self.atomic_numbers[n] == 1)
            .count()
    }

    fn degree(&self, i: usize) -> usize {
        self.adj[i].len()
    }

    fn valence(&self, i: usize) -> usize {
        self.bonds
            .iter()
            .enumerate()
            .filter(|(_, &(a, b))| a == i || b == i)
            .map(|(bi, _)| usize::from(self.bond_order[bi]))
            .sum()
    }
}

/// RDKit writes hybridization as `^1` = SP, `^2` = SP2, `^3` = SP3.
fn hyb_matches(h: Hybridization, n: u8) -> bool {
    matches!(
        (n, h),
        (1, Hybridization::Sp)
            | (2, Hybridization::Sp2)
            | (3, Hybridization::Sp3)
            | (4, Hybridization::Sp3d)
            | (5, Hybridization::Sp3d2)
    )
}

/// A molecule prepared for matching: its view, its graph built once, and a memo of each recursive
/// subquery's result.
///
/// The torsion library is matched pattern by pattern against the same molecule, and many patterns
/// carry `$(...)` recursive predicates that are evaluated at many candidate atoms. Rebuilding the
/// molecule graph per pattern, or re-running a recursive subquery per candidate atom, is what makes
/// the naive matcher quadratic; both are done once here instead. RDKit caches recursive query
/// results on the molecule in exactly the same spirit.
pub struct Matcher<'a> {
    pub view: MolView<'a>,
    graph: vf2::Graph,
    /// per recursive subpattern (keyed by its address), the atoms that can be its first atom
    rec_roots: std::cell::RefCell<HashMap<usize, std::rc::Rc<Vec<bool>>>>,
}

impl<'a> Matcher<'a> {
    pub fn new(view: MolView<'a>) -> Self {
        let graph = vf2::Graph::new(view.atomic_numbers.len(), view.bonds);
        Matcher {
            view,
            graph,
            rec_roots: std::cell::RefCell::new(HashMap::new()),
        }
    }

    /// Whether the recursive subpattern `sub` can root at atom `i`, computing its root set once.
    fn recursive_root(&self, sub: &Pattern, i: usize) -> bool {
        let key = sub as *const Pattern as usize;
        if let Some(bits) = self.rec_roots.borrow().get(&key) {
            return bits[i];
        }
        // Not cached: match the whole subpattern once and record every atom its first atom lands
        // on. No cache borrow is held across this, so nested recursive queries cache independently.
        let mut bits = vec![false; self.view.atomic_numbers.len()];
        for m in find_matches(sub, self) {
            bits[m[0]] = true;
        }
        let bits = std::rc::Rc::new(bits);
        let hit = bits[i];
        self.rec_roots.borrow_mut().insert(key, bits);
        hit
    }
}

/// Does atom `i` satisfy `pred`?
pub fn atom_matches(pred: &AtomPred, m: &Matcher, i: usize) -> bool {
    let mol = &m.view;
    match pred {
        AtomPred::Any => true,
        AtomPred::Symbol { z, aromatic } => {
            mol.atomic_numbers[i] == *z
                && aromatic.is_none_or(|want| mol.aromatic[i] == want)
        }
        AtomPred::AnyAromatic => mol.aromatic[i],
        AtomPred::AnyAliphatic => !mol.aromatic[i],
        AtomPred::AtomicNum(z) => mol.atomic_numbers[i] == *z,
        // `H` with no digit is exactly one hydrogen, not "at least one":
        // smarts.yy's H_TOKEN production builds makeAtomHCountQuery(1).
        AtomPred::TotalH(None) => mol.total_h(i) == 1,
        AtomPred::TotalH(Some(n)) => mol.total_h(i) == usize::from(*n),
        AtomPred::Connectivity(n) => mol.degree(i) == usize::from(*n),
        AtomPred::Degree(n) => mol.degree(i) == usize::from(*n),
        AtomPred::Valence(n) => mol.valence(i) == usize::from(*n),
        AtomPred::RingCount(None) => mol.ring_count[i] >= 1,
        AtomPred::RingCount(Some(n)) => mol.ring_count[i] == usize::from(*n),
        // RDKit's bare `r<n>` matches iff the atom's SMALLEST ring is exactly size n — not
        // membership in any ring of size n. A fusion atom shared between a five- and six-ring
        // (smallest ring 5) fails `r6` even though it is topologically in a six-ring; matching it
        // would over-select. Verified against the RDKit substructure matcher on indane and
        // hydrindane (fusion smallest-ring 5 → matches r5 only) and naphthalene (smallest-ring
        // 6 → matches r6). The Python `RingInfo` APIs (`IsAtomInRingOfSize`, `AtomRingSizes`)
        // report full membership, so they disagree with the matcher on fused atoms — the matcher
        // keys on the minimum, and it is the matcher the torsion library is applied through.
        AtomPred::RingSize(n) => mol.ring_sizes[i]
            .iter()
            .min()
            .is_some_and(|smallest| *smallest == usize::from(*n)),
        // RDKit's r{n-} asks whether the SMALLEST ring containing the atom is at least n — that
        // is, the atom is in no ring smaller than n. Not "is in some ring of at least n":
        // a macrocycle atom fused to a six-ring fails this, and matching it would over-select.
        // (QueryOps.h: queryAtomMinRingSize -> RingInfo::minAtomRingSize)
        AtomPred::MinRingSize(n) => mol
            .ring_sizes[i]
            .iter()
            .min()
            .is_some_and(|smallest| *smallest >= usize::from(*n)),
        AtomPred::RingBondCount(n) => mol.ring_bonds_at[i] == usize::from(*n),
        AtomPred::Hybridization(n) => hyb_matches(mol.hybridization[i], *n),
        AtomPred::Charge(c) => mol.charges[i] == *c,
        AtomPred::Recursive(sub) => m.recursive_root(sub, i),
        AtomPred::Not(p) => !atom_matches(p, m, i),
        AtomPred::And(ps) => ps.iter().all(|p| atom_matches(p, m, i)),
        AtomPred::Or(ps) => ps.iter().any(|p| atom_matches(p, m, i)),
    }
}

/// Does the bond at index `bi` satisfy `pred`?
pub fn bond_matches(pred: &BondPred, mol: &MolView, bi: usize) -> bool {
    match pred {
        // an unwritten bond is single or aromatic
        BondPred::Default => mol.bond_order[bi] == 1 || mol.bond_aromatic[bi],
        BondPred::Any => true,
        BondPred::Single => mol.bond_order[bi] == 1 && !mol.bond_aromatic[bi],
        BondPred::Double => mol.bond_order[bi] == 2 && !mol.bond_aromatic[bi],
        BondPred::Triple => mol.bond_order[bi] == 3,
        BondPred::Aromatic => mol.bond_aromatic[bi],
        BondPred::Ring => mol.bond_in_ring[bi],
        BondPred::Up | BondPred::Down => mol.bond_order[bi] == 1,
        BondPred::Not(p) => !bond_matches(p, mol, bi),
        BondPred::And(ps) => ps.iter().all(|p| bond_matches(p, mol, bi)),
        BondPred::Or(ps) => ps.iter().any(|p| bond_matches(p, mol, bi)),
    }
}

/// `SubstructMatchParams::maxMatches`, the cap `getExperimentalTorsions` matches under.
const MAX_MATCHES: usize = 1000;

/// Build a pattern's query graph. Molecule-independent, so it is computed once per pattern and
/// reused across every molecule.
pub fn pattern_graph(pattern: &Pattern) -> vf2::Graph {
    let edges: Vec<(usize, usize)> = pattern.bonds.iter().map(|&(a, b, _)| (a, b)).collect();
    vf2::Graph::new(pattern.atoms.len(), &edges)
}

/// All mappings of `pattern` onto the molecule, using a prebuilt query graph.
///
/// In discovery order and not deduplicated, because the caller keeps the first match for a given
/// central bond and discards the rest.
pub fn find_matches_with(pattern: &Pattern, g_pat: &vf2::Graph, m: &Matcher) -> Vec<Vec<usize>> {
    vf2::all_matches(
        g_pat,
        &m.graph,
        |pi, mi| atom_matches(&pattern.atoms[pi], m, mi),
        |pk, mb| bond_matches(&pattern.bonds[pk].2, &m.view, mb),
        MAX_MATCHES,
    )
}

/// All mappings of `pattern` onto the molecule, building the query graph on the fly.
///
/// Used for recursive subqueries and in tests; the hot per-molecule loop uses
/// [`find_matches_with`] with a graph built once.
pub fn find_matches(pattern: &Pattern, m: &Matcher) -> Vec<Vec<usize>> {
    find_matches_with(pattern, &pattern_graph(pattern), m)
}

/// Per-atom and per-bond ring data, derived once from a ring set.
pub fn ring_data(
    n_atoms: usize,
    bonds: &[(usize, usize)],
    rings: &[Vec<usize>],
    in_ring: &[bool],
) -> (Vec<Vec<usize>>, Vec<usize>, Vec<usize>) {
    let mut sizes = vec![Vec::new(); n_atoms];
    let mut count = vec![0usize; n_atoms];
    for r in rings {
        for &a in r {
            if a < n_atoms {
                sizes[a].push(r.len());
                count[a] += 1;
            }
        }
    }
    let mut ring_bonds_at = vec![0usize; n_atoms];
    for (bi, &(a, b)) in bonds.iter().enumerate() {
        if in_ring.get(bi).copied().unwrap_or(false) {
            ring_bonds_at[a] += 1;
            ring_bonds_at[b] += 1;
        }
    }
    (sizes, count, ring_bonds_at)
}

/// Assemble everything a match needs from a SMILES, for tests and for the torsion pass.
pub struct Perceived {
    pub atomic_numbers: Vec<u8>,
    pub charges: Vec<i8>,
    pub aromatic: Vec<bool>,
    pub hybridization: Vec<Hybridization>,
    pub bonds: Vec<(usize, usize)>,
    pub bond_order: Vec<u8>,
    pub bond_aromatic: Vec<bool>,
    pub bond_in_ring: Vec<bool>,
    pub rings: Vec<Vec<usize>>,
    /// per atom, `getValence(EXPLICIT) + getValence(IMPLICIT)` — RDKit's `Atom::getTotalValence`,
    /// which the UFF atom typer reads for oxidation-state suffixes. Invariant under `AddHs`.
    pub total_valence: Vec<i32>,
    /// per atom, whether any incident bond is conjugated — RDKit's `atomHasConjugatedBond`.
    pub atom_conjugated: Vec<bool>,
    /// per bond, whether it is conjugated — RDKit's `Bond::getIsConjugated`. Bonds to appended
    /// hydrogens are never conjugated.
    pub bond_conjugated: Vec<bool>,
    /// per bond, RDKit's `Bond::getStereo` (double-bond E/Z stereo).
    pub bond_stereo: Vec<Stereo>,
    /// per bond, RDKit's two `getStereoAtoms` reference atoms (`-1` when unset).
    pub stereo_atoms: Vec<[i64; 2]>,
    /// per atom, RDKit's tetrahedral `getChiralTag`.
    pub chiral_tags: Vec<crate::chirality::ChiralTag>,
}

/// Double-bond stereo, mirroring RDKit's `Bond::BondStereo`. `set14` treats `Z`≡`Cis` and
/// `E`≡`Trans`; the `> Any` threshold marks a definite assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stereo {
    None,
    Any,
    Z,
    E,
    Cis,
    Trans,
}

/// Perceive a SMILES all the way to the post-`AddHs` molecule the matcher sees.
pub fn perceive(smiles: &str) -> Result<Perceived, crate::ParseError> {
    use yowl::feature::BondKind;
    // Mirror RDKit's sanitizeMol order exactly, so everything downstream sees PERCEIVED aromaticity
    // (RDKit re-perceives regardless of how the SMILES wrote the rings): parse → valence#1 → SSSR →
    // Kekulize → assignRadicals → setAromaticity → setConjugation/Hybridization → valence#2.
    let g0 = crate::parse(smiles)?;
    let n_heavy = g0.atoms.len();

    // step 2 — first valence pass (as-written bonds); feeds kekulization
    let counts1 = crate::hydrogens::counts(&g0);
    // step 3 — SSSR
    let rings = crate::sssr::symmetrized_sssr_masked(n_heavy, &g0.bonds, &g0.dative_donor);

    // step 4 — Kekulize: as-written per-bond order/aromaticity in, definite single/double/triple out
    let atomic_heavy: Vec<u8> = g0.atoms.iter().map(|a| a.atomic_number).collect();
    let charges_heavy: Vec<i8> = g0.atoms.iter().map(|a| a.charge).collect();
    let arom_asw: Vec<bool> = g0.atoms.iter().map(|a| a.aromatic_as_written).collect();
    let tv1: Vec<i32> =
        (0..n_heavy).map(|i| counts1.explicit_valence[i] + counts1.implicit_valence[i]).collect();
    let mut bo_asw = vec![1u8; g0.bonds.len()];
    let mut ba_asw = vec![false; g0.bonds.len()];
    for (bi, &(a, b)) in g0.bonds.iter().enumerate() {
        let both = arom_asw[a] && arom_asw[b];
        let c = crate::hydrogens::contribution(g0.bond_kinds[bi], both);
        ba_asw[bi] = (c - 1.5).abs() < 1e-9;
        bo_asw[bi] = if (c - 2.0).abs() < 1e-9 { 2 } else if (c - 3.0).abs() < 1e-9 { 3 } else { 1 };
    }
    let kekule = crate::kekulize::kekule_bond_orders_raw(&crate::kekulize::KekMol {
        atomic_numbers: &atomic_heavy,
        bonds: &g0.bonds,
        bond_order: &bo_asw,
        bond_aromatic: &ba_asw,
        aromatic: &arom_asw,
        charges: &charges_heavy,
        total_valence: &tv1,
        total_num_hs: &counts1.total_num_hs,
    });

    // step 5 (assignRadicals — counts1.radicals is invariant vs the kekulized form) + step 6
    // (setAromaticity): perceive aromaticity from the kekulized structure
    let (atom_arom, bond_arom) = crate::aromaticity::perceive(
        &crate::aromaticity::AromInput {
            atomic_number: &atomic_heavy,
            charge: &charges_heavy,
            radicals: &counts1.radicals,
            total_num_hs: &counts1.total_num_hs,
            bonds: &g0.bonds,
            bond_orders: &kekule,
        },
        &rings,
    );

    // Rewrite the graph to RDKit's post-setAromaticity state: aromatic bonds → Aromatic, the rest →
    // their kekulized order. Then `counts` and `hybrid` (which read bond kinds + atom flags) see the
    // perceived state transparently — biphenyl's inter-ring single bond stays single since
    // `contribution` only treats an *Elided* bond as aromatic, never an explicit Single.
    let mut g = g0.clone();
    for (i, a) in g.atoms.iter_mut().enumerate() {
        a.aromatic_as_written = atom_arom[i];
    }
    for (bi, k) in g.bond_kinds.iter_mut().enumerate() {
        // Preserve directional single bonds (`/`,`\` = Up/Down): they are single bonds whose direction
        // encodes double-bond (E/Z) stereo, which is read separately downstream. RDKit keeps bond
        // direction distinct from bond type; rewriting these to plain Single would erase the stereo.
        if matches!(*k, BondKind::Up | BondKind::Down) {
            continue;
        }
        *k = if bond_arom[bi] {
            BondKind::Aromatic
        } else {
            match kekule[bi] {
                2 => BondKind::Double,
                3 => BondKind::Triple,
                _ => BondKind::Single,
            }
        };
    }

    // step 7 (conjugation + hybridization) + step 8 (second valence pass), all against perceived state
    let counts2 = crate::hydrogens::counts(&g);
    // adjustHs (RDKit `MolOps::adjustHs`): aromatization can drop an atom's implicit H — e.g. a
    // pyrrole-type N reaches valence 3 from two aromatic bonds and would lose its H — so RDKit adds
    // the lost implicit Hs back as EXPLICIT Hs, preserving the total. Total H = bracket + max(the
    // first-pass implicit, the perceived-pass implicit).
    let mut counts = counts2.clone();
    for i in 0..n_heavy {
        let bracket = counts1.total_num_hs[i] - counts1.implicit_valence[i];
        counts.total_num_hs[i] =
            bracket + counts1.implicit_valence[i].max(counts2.implicit_valence[i]);
    }
    let ex = crate::addhs::add_hs(&g, &counts);
    let (conj, hyb) = crate::hybrid::perceive(&g, &counts);

    let n = ex.atomic_numbers.len();
    let mut charges = vec![0i8; n];
    let mut aromatic = vec![false; n];
    let mut hybridization = vec![Hybridization::Unspecified; n];
    for (i, a) in g.atoms.iter().enumerate() {
        charges[i] = a.charge;
        aromatic[i] = a.aromatic_as_written;
        hybridization[i] = hyb[i];
    }

    let mut bond_order = Vec::with_capacity(ex.bonds.len());
    let mut bond_aromatic = Vec::with_capacity(ex.bonds.len());
    for (bi, &(a, b)) in ex.bonds.iter().enumerate() {
        if bi < g.bond_kinds.len() {
            let ar = aromatic[a] && aromatic[b];
            let c = crate::hydrogens::contribution(g.bond_kinds[bi], ar);
            // Bond aromaticity is the perceived value (RDKit's markAtomsBondsArom), not re-derived
            // from the bond kind: an aromatic bond that is also a stereo-direction bond (`/`,`\`)
            // keeps its Up/Down kind — so `contribution` would read it as single — but it is still
            // aromatic. RDKit stores bond type (AROMATIC) and direction as independent fields.
            bond_aromatic.push(bond_arom[bi]);
            bond_order.push(if (c - 2.0).abs() < 1e-9 {
                2
            } else if (c - 3.0).abs() < 1e-9 {
                3
            } else {
                1
            });
        } else {
            // an appended hydrogen bond
            bond_aromatic.push(false);
            bond_order.push(1);
        }
    }
    let n_bonds = ex.bonds.len();
    let heavy_in_ring = crate::sssr::ring_bond_flags_masked(g.atoms.len(), &g.bonds, &g.dative_donor);
    let mut bond_in_ring = heavy_in_ring;
    bond_in_ring.resize(n_bonds, false);

    // getTotalValence = explicit + implicit valence, invariant under AddHs. Heavy atoms keep their
    // pre-AddHs index, so the counts line up; each appended hydrogen has total valence 1.
    let n_heavy = g.atoms.len();
    let mut total_valence = vec![1i32; n];
    for (i, tv) in total_valence.iter_mut().take(n_heavy).enumerate() {
        *tv = counts.explicit_valence[i] + counts.implicit_valence[i];
    }
    // atomHasConjugatedBond: an atom is conjugated if any incident (heavy) bond is. Bonds to
    // appended hydrogens are never conjugated, so only the heavy bonds contribute.
    let mut atom_conjugated = vec![false; n];
    let mut bond_conjugated = vec![false; ex.bonds.len()];
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        if conj[bi] {
            atom_conjugated[a] = true;
            atom_conjugated[b] = true;
            bond_conjugated[bi] = true;
        }
    }

    let mut out = Perceived {
        atomic_numbers: ex.atomic_numbers,
        charges,
        aromatic,
        hybridization,
        bonds: ex.bonds,
        bond_order,
        bond_aromatic,
        bond_in_ring,
        rings,
        total_valence,
        atom_conjugated,
        bond_conjugated,
        // double-bond stereo is filled by stereo::assign below
        bond_stereo: vec![Stereo::None; n_bonds],
        stereo_atoms: vec![[-1, -1]; n_bonds],
        chiral_tags: {
            // tetrahedral tags are perceived on the heavy graph; heavy atoms keep their index post
            // AddHs, and appended hydrogens are never stereocenters
            let heavy_bond_orders: Vec<f64> = g
                .bonds
                .iter()
                .enumerate()
                .map(|(bi, &(a, b))| {
                    let both = g.atoms[a].aromatic_as_written && g.atoms[b].aromatic_as_written;
                    crate::hydrogens::contribution(g.bond_kinds[bi], both)
                })
                .collect();
            let heavy = crate::chirality::compute_tags(
                &g.chirality,
                &g.bonds,
                &heavy_bond_orders,
                &counts.implicit_valence,
            );
            let mut tags = vec![crate::chirality::ChiralTag::None; n];
            tags[..heavy.len()].copy_from_slice(&heavy);
            tags
        },
    };
    crate::stereo::assign(&g, &mut out);
    Ok(out)
}

impl Perceived {
    /// Ring membership over bonds: for each ring of `rings`, in the same order, the indices of the ring
    /// bonds whose two atoms both lie in it.
    pub fn bond_rings(&self) -> Vec<Vec<usize>> {
        self.rings
            .iter()
            .map(|ring| {
                let set: std::collections::HashSet<usize> = ring.iter().copied().collect();
                self.bonds
                    .iter()
                    .enumerate()
                    .filter(|&(bi, &(a, b))| {
                        set.contains(&a) && set.contains(&b) && self.bond_in_ring.get(bi).copied().unwrap_or(false)
                    })
                    .map(|(bi, _)| bi)
                    .collect()
            })
            .collect()
    }

    pub fn view(&self) -> MolView<'_> {
        let n = self.atomic_numbers.len();
        let mut adj = vec![Vec::new(); n];
        for &(a, b) in &self.bonds {
            adj[a].push(b);
            adj[b].push(a);
        }
        let (ring_sizes, ring_count, ring_bonds_at) =
            ring_data(n, &self.bonds, &self.rings, &self.bond_in_ring);
        MolView {
            atomic_numbers: &self.atomic_numbers,
            charges: &self.charges,
            aromatic: &self.aromatic,
            hybridization: &self.hybridization,
            adj,
            bonds: &self.bonds,
            bond_order: &self.bond_order,
            bond_aromatic: &self.bond_aromatic,
            bond_in_ring: &self.bond_in_ring,
            ring_sizes,
            ring_count,
            ring_bonds_at,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// cyclo(L-Ala)6: six macrocyclic amides in an 18-membered ring.
    const CYCLO_ALA6: &str = "C[C@@H]1C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N1";

    #[test]
    fn finds_a_simple_pattern() {
        let p = perceive("CCO").expect("perceives");
        let m = Matcher::new(p.view());
        let pat = crate::smarts::parse("[CX4:1][OX2:2]").expect("parses");
        assert_eq!(find_matches(&pat, &m).len(), 1, "one C-O match expected");
    }

    #[test]
    fn ring_size_predicate_selects_the_macrocycle() {
        let p = perceive(CYCLO_ALA6).expect("perceives");
        let mt = Matcher::new(p.view());
        // any carbon in a ring of at least 9
        let big = crate::smarts::parse("[C;r{9-}:1]").expect("parses");
        let small = crate::smarts::parse("[C;r{19-}:1]").expect("parses");
        let n_big = find_matches(&big, &mt).len();
        assert!(n_big > 0, "the 18-ring should satisfy r{{9-}}");
        assert_eq!(find_matches(&small, &mt).len(), 0, "nothing is in a 19+ ring");
    }

    #[test]
    fn bare_ring_size_keys_on_the_smallest_ring() {
        // Indane: two carbons are shared between the five- and six-membered rings. RDKit's bare
        // `r<n>` keys on the smallest ring, so those fusion atoms match `r5` but NOT `r6`, even
        // though they are topologically in the six-ring. A membership reading would wrongly match
        // `r6` and pick the wrong torsion pattern for aryl-carbonyl bonds off fused heteroaromatics.
        let p = perceive("C1CCc2ccccc21").expect("perceives");
        let mt = Matcher::new(p.view());
        let r5 = find_matches(&crate::smarts::parse("[*;r5:1]").expect("parses"), &mt).len();
        let r6 = find_matches(&crate::smarts::parse("[*;r6:1]").expect("parses"), &mt).len();
        // 5 atoms have smallest-ring 5 (the two fusion + three CH2); 4 have smallest-ring 6.
        assert_eq!(r5, 5, "five atoms have a smallest ring of size 5");
        assert_eq!(r6, 4, "only the four non-fusion aromatic carbons have smallest ring 6");
    }

    #[test]
    fn divyas_patched_amide_pattern_matches() {
        let p = perceive(CYCLO_ALA6).expect("perceives");
        let mt = Matcher::new(p.view());
        // the patched pattern: a macrocyclic amide with force constant 80.0
        let pat = crate::smarts::parse("[C:1][C;r{9-}:2](=O)@;-[NX3H1;r:3][CX4H1:4]")
            .expect("parses");
        let m = find_matches(&pat, &mt);
        assert!(!m.is_empty(), "cyclo(L-Ala)6 must match the patched amide torsion");
        println!("  patched amide pattern matched {} times", m.len());
    }
}
