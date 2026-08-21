//! The SMARTS abstract syntax tree, and lowering from the parse tree.
//!
//! Predicates are kept as a tree rather than flattened, so `!`, `&`, `,` and `;` retain their
//! precedence and evaluation is a straight recursive walk. Every primitive here is one the torsion
//! library actually uses; the grammar rejects anything else.

/// A condition on one atom.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomPred {
    /// `*`
    Any,
    /// An element symbol; `aromatic` distinguishes `c` from `C`.
    Symbol { z: u8, aromatic: Option<bool> },
    /// `a`
    AnyAromatic,
    /// `A`
    AnyAliphatic,
    /// `#N`
    AtomicNum(u8),
    /// `H` (at least one) or `H<n>` (exactly n)
    TotalH(Option<u8>),
    /// `X<n>` — total connections including hydrogens
    Connectivity(u8),
    /// `D<n>` — explicit connections
    Degree(u8),
    /// `v<n>` — total bond order
    Valence(u8),
    /// `R` (in any ring) or `R<n>` (in exactly n rings)
    RingCount(Option<u8>),
    /// `r<n>` — in a ring of exactly this size
    RingSize(u8),
    /// `r{n-}` — in a ring of at least this size
    MinRingSize(u8),
    /// `x<n>` — number of ring bonds
    RingBondCount(u8),
    /// `^n` — hybridization
    Hybridization(u8),
    Charge(i8),
    /// `$(...)` — matches when the sub-pattern matches with its first atom here
    Recursive(Box<Pattern>),
    Not(Box<AtomPred>),
    And(Vec<AtomPred>),
    Or(Vec<AtomPred>),
}

/// A condition on one bond.
#[derive(Debug, Clone, PartialEq)]
pub enum BondPred {
    /// No bond symbol written: single or aromatic.
    Default,
    /// `~`
    Any,
    Single,
    Double,
    Triple,
    Aromatic,
    /// `@` — a ring bond
    Ring,
    Up,
    Down,
    Not(Box<BondPred>),
    And(Vec<BondPred>),
    Or(Vec<BondPred>),
}

/// A parsed pattern: atoms in parse order, their map numbers, and the bonds between them.
#[derive(Debug, Clone, PartialEq)]
pub struct Pattern {
    pub atoms: Vec<AtomPred>,
    pub maps: Vec<Option<u8>>,
    pub bonds: Vec<(usize, usize, BondPred)>,
}

impl Pattern {
    /// Atom index carrying a given map number, if any.
    pub fn atom_with_map(&self, map: u8) -> Option<usize> {
        self.maps.iter().position(|m| *m == Some(map))
    }

    /// The torsion quadruple, by map numbers 1..4.
    pub fn torsion_quad(&self) -> Option<[usize; 4]> {
        Some([
            self.atom_with_map(1)?,
            self.atom_with_map(2)?,
            self.atom_with_map(3)?,
            self.atom_with_map(4)?,
        ])
    }
}
