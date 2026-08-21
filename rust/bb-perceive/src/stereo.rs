//! Double-bond stereo from SMILES directional bonds.
//!
//! A port of RDKit's `MolOps::setBondStereoFromDirections` (`Chirality.cpp`), which `SmilesToMol`
//! calls to turn the `/` and `\` markers into `Bond::getStereo`. For each double bond it finds a
//! directed single bond at each end, records their far atoms as the stereo reference atoms, and sets
//! CIS or TRANS from whether the two (orientation-normalised) directions agree.
//!
//! The directional markers come straight from yowl (`BondKind::Up` = `/`, `Down` = `\`), and our
//! bond list carries RDKit's begin→end orientation (verified by the AddHs bond-sequence gate), so
//! the begin/end tests below line up with RDKit's. `set14` treats CIS≡Z and TRANS≡E, so producing
//! the CIS/TRANS representation rather than RDKit's modern Z/E is equivalent for the bounds matrix.

use crate::smarts_match::{Perceived, Stereo};
use crate::SmilesGraph;
use yowl::feature::BondKind;

/// The direction of a bond as written, `true` for `/` (up), `false` for `\` (down).
fn direction(kind: BondKind) -> Option<bool> {
    match kind {
        BondKind::Up => Some(true),
        BondKind::Down => Some(false),
        _ => None,
    }
}

/// Assign double-bond stereo into `out`, from the directional bonds of the heavy-atom graph `g`.
///
/// Heavy bonds occupy the first `g.bonds.len()` entries of the post-`AddHs` bond list in the same
/// order, so bond indices carry straight over.
pub fn assign(g: &SmilesGraph, out: &mut Perceived) {
    let n = g.atoms.len();
    // incident heavy bonds per atom, in bond-index order (matching RDKit's getAtomBonds order)
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        adj[a].push(bi);
        adj[b].push(bi);
    }

    // the first non-double directed bond at an atom — RDKit's getNeighboringDirectedBond
    let neighboring_directed = |atom: usize| -> Option<usize> {
        adj[atom]
            .iter()
            .copied()
            .find(|&bi| out.bond_order[bi] != 2 && direction(g.bond_kinds[bi]).is_some())
    };

    for (bi, &(begin, end)) in g.bonds.iter().enumerate() {
        if out.bond_order[bi] != 2 || out.bond_stereo[bi] == Stereo::Any {
            continue; // only definite (non-Any) double bonds
        }
        let (Some(db), Some(de)) = (neighboring_directed(begin), neighboring_directed(end)) else {
            continue; // needs a directed bond at each end
        };

        let (dbb, dbe) = g.bonds[db]; // directed bond at the begin end
        let (deb, dee) = g.bonds[de]; // directed bond at the end end
        let begin_side_atom = if dbb == begin { dbe } else { dbb };
        let end_side_atom = if deb == end { dee } else { deb };
        out.stereo_atoms[bi] = [begin_side_atom as i64, end_side_atom as i64];

        // normalise each direction to point away from the double-bond atom, exactly as RDKit does:
        // flip the begin-side dir if that bond *begins* at the stereo atom, the end-side dir if it
        // *ends* at the stereo atom.
        let mut begin_dir = direction(g.bond_kinds[db]).expect("directed");
        if dbb == begin {
            begin_dir = !begin_dir;
        }
        let mut end_dir = direction(g.bond_kinds[de]).expect("directed");
        if dee == end {
            end_dir = !end_dir;
        }
        out.bond_stereo[bi] = if begin_dir == end_dir { Stereo::Trans } else { Stereo::Cis };
    }
}
