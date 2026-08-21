//! Tetrahedral chirality perception: the per-atom chiral tag RDKit's SMILES parser assigns.
//!
//! Reproduces `SmilesParseOps::AdjustAtomChiralityFlags`. The raw `@`/`@@` is read from a dedicated
//! [`Follower`] (which, unlike yowl's `Builder`, receives the configuration *before* yowl's internal
//! `invert_configuration`), then adjusted to the built-molecule bond order:
//!
//! 1. `GetBondOrdering` — the SMILES order of the atom's heavy bonds: non-closure neighbours sorted
//!    by atom index, with the ring-closure bonds inserted at the atom's own position (i.e. after any
//!    lower-indexed neighbour, before the higher-indexed ones), in ring-open order.
//! 2. `getPerturbationOrder` — the swap parity between that SMILES order and the built bond order
//!    (`countSwapsToInterconvert`).
//! 3. `chiralAtomNeedsTagInversion` — a degree-3 correction (implicit-H-first, single ring closure).
//!
//! An odd total swap count inverts the tag. `@` (yowl `TH1`) starts as CCW, `@@` (`TH2`) as CW.

use yowl::feature::{AtomKind, BondKind, Configuration, Rnum};
use yowl::walk::Follower;

/// RDKit's `Atom::ChiralType` for the tetrahedral cases (`CHI_UNSPECIFIED=0`, `CW=1`, `CCW=2`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChiralTag {
    None,
    Cw,
    Ccw,
}

impl ChiralTag {
    fn inverted(self) -> Self {
        match self {
            ChiralTag::Cw => ChiralTag::Ccw,
            ChiralTag::Ccw => ChiralTag::Cw,
            ChiralTag::None => ChiralTag::None,
        }
    }
}

/// Per-atom parse data needed for chirality, captured in walk (atom-index) order — the same order
/// yowl's `Builder` and this crate's other followers number atoms in.
#[derive(Default, Clone, Debug, PartialEq)]
pub struct ChiralityData {
    next_id: usize,
    stack: Vec<usize>,
    /// open ring bonds: `rnum -> (opener atom, slot in that atom's closure list)`
    opens: std::collections::HashMap<Rnum, (usize, usize)>,
    /// raw `@`/`@@` per atom (pre yowl inversion)
    config: Vec<Option<Configuration>>,
    /// whether the atom began a fragment (RDKit's `_SmilesStart`)
    is_root: Vec<bool>,
    /// bracket hydrogen count per atom
    num_hs: Vec<u8>,
    /// ring-closure partner atoms per atom, in ring-open (written) order (`usize::MAX` until closed)
    closure_partners: Vec<Vec<usize>>,
}

const UNRESOLVED: usize = usize::MAX;

fn config_and_hs(kind: &AtomKind) -> (Option<Configuration>, u8) {
    match kind {
        AtomKind::Bracket { configuration, hcount, .. } => {
            (*configuration, hcount.as_ref().map_or(0, u8::from))
        }
        AtomKind::Symbol(_) => (None, 0),
    }
}

impl ChiralityData {
    /// Whether atom `i` carries a parsed `@`/`@@` tag.
    pub fn is_chiral(&self, i: usize) -> bool {
        self.config[i].is_some()
    }

    /// `removeHs`: one more explicit H folded onto atom `i`.
    pub fn bump_num_hs(&mut self, i: usize) {
        self.num_hs[i] += 1;
    }

    /// `removeHs` odd perturbation: swap the raw `@`↔`@@` sense on atom `i`.
    pub fn invert(&mut self, i: usize) {
        self.config[i] = self.config[i].map(|c| match c {
            Configuration::TH1 => Configuration::TH2,
            Configuration::TH2 => Configuration::TH1,
            other => other,
        });
    }

    /// Compact the per-atom vectors to the surviving atoms (`new_index[old] = new`, `usize::MAX` for
    /// removed) after `removeHs`. Closure-partner atom indices are remapped too; a partner pointing
    /// at a removed atom cannot occur (a directional/closure bond keeps its H).
    pub fn compact(&mut self, new_index: &[usize], new_len: usize) {
        let mut config = vec![None; new_len];
        let mut is_root = vec![false; new_len];
        let mut num_hs = vec![0u8; new_len];
        let mut closure_partners: Vec<Vec<usize>> = vec![Vec::new(); new_len];
        for old in 0..new_index.len() {
            let ni = new_index[old];
            if ni == usize::MAX {
                continue;
            }
            config[ni] = self.config[old];
            is_root[ni] = self.is_root[old];
            num_hs[ni] = self.num_hs[old];
            closure_partners[ni] = self.closure_partners[old]
                .iter()
                .map(|&p| if p == UNRESOLVED { UNRESOLVED } else { new_index[p] })
                .collect();
        }
        self.config = config;
        self.is_root = is_root;
        self.num_hs = num_hs;
        self.closure_partners = closure_partners;
    }

    fn add_atom(&mut self, kind: &AtomKind, is_root: bool) -> usize {
        let (cfg, hs) = config_and_hs(kind);
        let id = self.next_id;
        self.next_id += 1;
        self.config.push(cfg);
        self.is_root.push(is_root);
        self.num_hs.push(hs);
        self.closure_partners.push(Vec::new());
        id
    }
}

impl Follower for ChiralityData {
    fn root(&mut self, kind: AtomKind) {
        let id = self.add_atom(&kind, true);
        self.stack.push(id);
    }

    fn extend(&mut self, _bond_kind: BondKind, kind: AtomKind) {
        let id = self.add_atom(&kind, false);
        self.stack.push(id);
    }

    fn join(&mut self, _bond_kind: BondKind, rnum: Rnum) {
        let sid = *self.stack.last().expect("stack is non-empty");
        match self.opens.remove(&rnum) {
            Some((opener, slot)) => {
                // closing: resolve the opener's placeholder and record the partner here
                self.closure_partners[opener][slot] = sid;
                self.closure_partners[sid].push(opener);
            }
            None => {
                let slot = self.closure_partners[sid].len();
                self.closure_partners[sid].push(UNRESOLVED);
                self.opens.insert(rnum, (sid, slot));
            }
        }
    }

    fn pop(&mut self, depth: usize) {
        for _ in 0..depth {
            self.stack.pop();
        }
    }
}

/// `countSwapsToInterconvert`: swaps to turn `probe` into `ref` (both permutations of one set).
fn count_swaps(reference: &[usize], probe: &[usize]) -> usize {
    let mut probe = probe.to_vec();
    let mut swaps = 0;
    for i in 0..reference.len() {
        if probe[i] != reference[i] {
            let j = (i + 1..probe.len())
                .find(|&j| probe[j] == reference[i])
                .expect("permutations of the same set");
            probe.swap(i, j);
            swaps += 1;
        }
    }
    swaps
}

/// Compute the per-atom chiral tags. `bonds` is the built heavy-bond list in final (bond-index)
/// order; `bond_orders` and `implicit_valence` describe each bond and atom for the inversion rule.
pub fn compute_tags(
    data: &ChiralityData,
    bonds: &[(usize, usize)],
    bond_orders: &[f64],
    implicit_valence: &[i32],
) -> Vec<ChiralTag> {
    let n = data.config.len();
    // built neighbour order: neighbours in the order their bonds appear in the final bond list
    let mut built: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut incident_unsaturated = vec![false; n];
    for (bi, &(a, b)) in bonds.iter().enumerate() {
        built[a].push(b);
        built[b].push(a);
        if bond_orders[bi] > 1.0 {
            incident_unsaturated[a] = true;
            incident_unsaturated[b] = true;
        }
    }

    (0..n)
        .map(|x| {
            let initial = match data.config[x] {
                Some(Configuration::TH1) => ChiralTag::Ccw,
                Some(Configuration::TH2) => ChiralTag::Cw,
                _ => return ChiralTag::None,
            };
            let closures = &data.closure_partners[x];
            // non-closure neighbours = built neighbours minus the closure partners (one bond each)
            let mut remaining = closures.clone();
            let mut non_closure: Vec<usize> = Vec::new();
            for &nbr in &built[x] {
                if let Some(pos) = remaining.iter().position(|&c| c == nbr) {
                    remaining.swap_remove(pos);
                } else {
                    non_closure.push(nbr);
                }
            }
            // GetBondOrdering: sorted lower-index neighbours, then closures (open order), then
            // sorted higher-index neighbours
            let mut lower: Vec<usize> = non_closure.iter().copied().filter(|&a| a < x).collect();
            let mut higher: Vec<usize> = non_closure.iter().copied().filter(|&a| a > x).collect();
            lower.sort_unstable();
            higher.sort_unstable();
            let mut smiles_order = lower;
            smiles_order.extend(closures.iter().copied());
            smiles_order.extend(higher);

            let mut nswaps = count_swaps(&built[x], &smiles_order);

            // chiralAtomNeedsTagInversion (degree-3)
            let degree = built[x].len();
            if degree == 3 {
                let num_hs = data.num_hs[x];
                let has_fourth_valence = num_hs == 1 || implicit_valence[x] == 1;
                if (data.is_root[x] && num_hs == 1)
                    || (!has_fourth_valence && closures.len() == 1 && !incident_unsaturated[x])
                {
                    nswaps += 1;
                }
            }

            if nswaps % 2 == 1 {
                initial.inverted()
            } else {
                initial
            }
        })
        .collect()
}
