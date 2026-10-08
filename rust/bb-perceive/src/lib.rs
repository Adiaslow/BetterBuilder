
//! Molecular perception in Rust, reproducing RDKit's: SMILES parsing, rings ([`rings`], [`sssr`]),
//! aromaticity, valence and hydrogens ([`valence`], [`hydrogens`], [`addhs`]), hybridization,
//! chirality and double-bond stereo, SMARTS matching ([`smarts`], [`smarts_match`], [`vf2`]), and the
//! experimental-torsion library ([`torsions`]) — each gated against RDKit.
//!
//! The crate root holds the parse: SMILES to a heavy-atom graph via yowl. It carries only what the
//! SMILES states — element, charge, bracket hydrogen count, bond order,
//! and aromaticity as written. Perceived aromaticity, implicit valence and hydrogen addition follow
//! RDKit's own algorithms and are not part of the parse.

pub mod cleanup;
pub mod recipe;
pub mod remove_hs;
pub mod rings;
pub mod smarts;
pub mod smarts_ast;
pub mod smarts_match;
pub mod chirality;
pub mod stereo;
pub mod sssr;
pub mod addhs;
pub mod angles;
pub mod aromaticity;
pub mod hybrid;
pub mod hydrogens;
pub mod kekulize;
pub mod torsion_lib;
pub mod torsions;
pub mod valence;
pub mod vf2;

use yowl::feature::{AtomKind, BondKind, Rnum, Symbol};
use yowl::graph::Builder;
use yowl::walk::Follower;

/// Records bond creation order the way RDKit's SMILES parser does.
///
/// `SmilesParse.cpp:256` calls `CloseMolRings` only after the whole string is parsed, so every
/// chain bond precedes every ring-closure bond; `CloseMolRings` then walks the atom bookmarks,
/// which are keyed by ring-closure digit, so closures are ordered by digit rather than by position
/// in the string. Atom numbering mirrors `yowl::graph::Builder` so the two agree on indices.
#[derive(Default)]
struct BondOrder {
    next_id: usize,
    stack: Vec<usize>,
    opens: std::collections::HashMap<Rnum, usize>,
    chain: Vec<(usize, usize)>,
    closures: Vec<(Rnum, usize, usize)>,
}

impl Follower for BondOrder {
    fn root(&mut self, _kind: AtomKind) {
        self.stack.push(self.next_id);
        self.next_id += 1;
    }

    fn extend(&mut self, _bond_kind: BondKind, _atom_kind: AtomKind) {
        let sid = *self.stack.last().expect("stack is non-empty");
        let tid = self.next_id;
        self.next_id += 1;
        self.chain.push((sid, tid));
        self.stack.push(tid);
    }

    fn join(&mut self, _bond_kind: BondKind, rnum: Rnum) {
        let sid = *self.stack.last().expect("stack is non-empty");
        match self.opens.remove(&rnum) {
            Some(tid) => self.closures.push((rnum, sid, tid)),
            None => {
                self.opens.insert(rnum, sid);
            }
        }
    }

    fn pop(&mut self, depth: usize) {
        for _ in 0..depth {
            self.stack.pop();
        }
    }
}


/// One heavy atom as the SMILES states it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeavyAtom {
    pub atomic_number: u8,
    pub charge: i8,
    /// Hydrogens written inside brackets; `None` for organic-subset atoms, whose count is implicit.
    pub bracket_hydrogens: Option<u8>,
    /// Aromatic as written in the SMILES, not as perceived.
    pub aromatic_as_written: bool,
    /// Whether the atom carries an explicit isotope (`[2H]`, `[13C]`, …). RDKit's `removeHs` keeps
    /// isotopic hydrogens (`removeIsotopes=false`), so deuterium/tritium must not be folded away.
    pub isotope: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SmilesGraph {
    pub atoms: Vec<HeavyAtom>,
    /// Sorted, deduplicated `(i, j)` with `i < j`.
    pub bonds: Vec<(usize, usize)>,
    pub bond_kinds: Vec<BondKind>,
    /// Per bond, the donor atom when `cleanUpOrganometallics` has made it a ligand→metal dative bond
    /// (else `None`). A dative bond is otherwise an ordinary single bond — `getBondTypeAsDouble` is
    /// 1.0, it counts toward degree and bond length as a single — but it contributes 0 to the
    /// donor's valence and electron count (RDKit's `Bond::DATIVE` / `getValenceContrib`).
    pub dative_donor: Vec<Option<usize>>,
    /// Raw tetrahedral-chirality parse data, for [`chirality::compute_tags`].
    pub chirality: chirality::ChiralityData,
}

#[derive(Debug)]
pub enum ParseError {
    Read(String),
    Graph(String),
    UnknownElement(String),
    /// An atom's explicit valence exceeds what its element allows — RDKit rejects such molecules
    /// during sanitization (`MolFromSmiles` returns null). Carries the offending atom index.
    Valence(usize),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::Read(s) => write!(f, "SMILES read error: {s}"),
            ParseError::Graph(s) => write!(f, "SMILES graph error: {s}"),
            ParseError::UnknownElement(s) => write!(f, "unknown element {s}"),
            ParseError::Valence(i) => write!(f, "atom {i} exceeds its allowed valence"),
        }
    }
}

impl std::error::Error for ParseError {}

/// yowl re-exports `mendeleev::Element`, so the atomic number comes straight from the periodic
/// table with no symbol round-trip. `*` (Star) has no element.
fn atomic_number(sym: &Symbol) -> Result<u8, ParseError> {
    match sym {
        Symbol::Aliphatic(e) | Symbol::Aromatic(e) => Ok(e.atomic_number() as u8),
        // `*` is a dummy atom (RDKit atomic number 0): no valence constraint, unspecified
        // hybridization, never an aromaticity candidate.
        Symbol::Star => Ok(0),
    }
}

/// Parse a SMILES into its heavy-atom graph.
pub fn parse(smiles: &str) -> Result<SmilesGraph, ParseError> {
    let mut builder = Builder::default();
    yowl::read::read(smiles, &mut builder, None).map_err(|e| ParseError::Read(format!("{e:?}")))?;
    let mut bond_order = BondOrder::default();
    yowl::read::read(smiles, &mut bond_order, None)
        .map_err(|e| ParseError::Read(format!("{e:?}")))?;
    let mut chirality = chirality::ChiralityData::default();
    yowl::read::read(smiles, &mut chirality, None)
        .map_err(|e| ParseError::Read(format!("{e:?}")))?;
    let graph = builder.build().map_err(|e| ParseError::Graph(format!("{e:?}")))?;

    let mut atoms = Vec::with_capacity(graph.len());
    for a in &graph {
        let (sym, hcount, charge, isotope) = match &a.kind {
            AtomKind::Symbol(s) => (s, None, 0i8, false),
            AtomKind::Bracket {
                symbol,
                hcount,
                charge,
                isotope,
                ..
            } => (
                symbol,
                // A bracket atom states its hydrogens explicitly (`noImplicit`); an absent count means
                // zero, NOT "implicit". `[O]` is atomic-oxygen-with-0-H (a radical), distinct from
                // organic `O` (water) — so brackets always yield `Some`, never `None`.
                Some(hcount.as_ref().map_or(0, u8::from)),
                charge.map_or(0, i8::from),
                isotope.is_some(),
            ),
        };
        atoms.push(HeavyAtom {
            atomic_number: atomic_number(sym)?,
            charge,
            bracket_hydrogens: hcount,
            aromatic_as_written: a.is_aromatic(),
            isotope,
        });
    }

    // Bond order is load-bearing, not cosmetic: ring perception walks an atom's neighbours in bond
    // order, and among equal-size rings that decides which one is found. Keep first-encounter order
    // rather than sorting.
    let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    let mut bonds: Vec<(usize, usize)> = Vec::new();
    let mut bond_kinds: Vec<BondKind> = Vec::new();
    for (i, a) in graph.iter().enumerate() {
        for b in &a.bonds {
            let (lo, hi) = if i < b.tid { (i, b.tid) } else { (b.tid, i) };
            if seen.insert((lo, hi)) {
                bonds.push((lo, hi));
                bond_kinds.push(b.kind);
            }
        }
    }

    // Bond order is load-bearing: ring perception walks an atom's neighbours in bond order, and
    // among equal-size rings that decides which is found. Rebuild the list in RDKit's order.
    let kinds: std::collections::HashMap<(usize, usize), BondKind> =
        bonds.iter().copied().zip(bond_kinds.iter().copied()).collect();
    let norm = |a: usize, b: usize| if a < b { (a, b) } else { (b, a) };
    let mut ordered: Vec<(usize, usize)> = bond_order.chain.clone();
    let mut closures = bond_order.closures.clone();
    closures.sort_by_key(|&(r, _, _)| r);
    ordered.extend(closures.iter().map(|&(_, a, b)| (a, b)));
    let bonds: Vec<(usize, usize)> = ordered;
    let bond_kinds: Vec<BondKind> = bonds
        .iter()
        .map(|&(a, b)| {
            let k = kinds.get(&norm(a, b)).copied().unwrap_or(BondKind::Elided);
            // `kinds` stores a bond's direction relative to `min→max` (first-seen from the lower
            // atom). Chain bonds are stored parent→child (min→max), but a ring closure is stored
            // closer→opener, i.e. reversed. `/` and `\` are direction-relative, so a reversed bond's
            // marker must flip — otherwise the E/Z sign of a double bond closed with a directional
            // ring bond (e.g. `.../2=N/O`) comes out inverted.
            if a > b {
                match k {
                    BondKind::Up => BondKind::Down,
                    BondKind::Down => BondKind::Up,
                    other => other,
                }
            } else {
                k
            }
        })
        .collect();

    let dative_donor = vec![None; bonds.len()];
    let graph = SmilesGraph {
        atoms,
        bonds,
        bond_kinds,
        dative_donor,
        chirality,
    };
    // MolFromSmiles removes explicit hydrogens (folding them into their neighbour's H count, fixing
    // stereocentres) before sanitizing, so they are re-added by AddHs at canonical indices.
    let mut graph = remove_hs::remove_hs(graph, bond_order.chain.len());
    // sanitizeMol step 1: cleanUp charge-separates nitro/azide/phosphoryl/hypervalent-halogen forms
    // before the valence check sees them (e.g. `OCl(=O)(=O)=O` → `[Cl+3]([O-])([O-])([O-])O`).
    cleanup::clean_up(&mut graph);
    // sanitizeMol step 2: cleanUpOrganometallics retypes ligand→metal single bonds from hypervalent
    // non-metals as dative, so the coordinated atom's valence stops counting the metal bond.
    cleanup::clean_up_organometallics(&mut graph);
    // Sanitization valence check: refuse molecules RDKit refuses (e.g. `[CH5]`, `[ClH2]`).
    if let Err(i) = hydrogens::check_valences(&graph) {
        return Err(ParseError::Valence(i));
    }
    Ok(graph)
}
