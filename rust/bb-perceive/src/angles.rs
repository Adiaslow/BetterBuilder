//! Bond-angle collection: the 1-3 layer of the distance bounds.
//!
//! Reproduces `DGeomHelpers::collectBondsAndAngles`. Bond pairs are visited as `(i, j)` with
//! `j > i` over bond index order, so the angle sequence follows bond order — which is why bond
//! order had to be right first. The central atom is whichever endpoint the two bonds share, and
//! the outer atoms keep the begin/end orientation the bonds were stored with.

use crate::addhs::Explicit;

/// An angle: outer, centre, outer, and whether it should be treated as linear.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Angle {
    pub atoms: [usize; 3],
    /// Set for a triple bond, or for two double bonds meeting at a two-coordinate centre —
    /// both cases where the three atoms are collinear.
    pub triple: bool,
}

/// Bond orders as the angle flag needs them, indexed like `bonds`.
pub fn collect(ex: &Explicit, bond_orders: &[u8], degrees: &[usize]) -> Vec<Angle> {
    let bonds = &ex.bonds;
    let mut out = Vec::new();
    for i in 0..bonds.len() {
        let (a11, a12) = bonds[i];
        for j in (i + 1)..bonds.len() {
            let (a21, a22) = bonds[j];
            if a11 != a21 && a11 != a22 && a12 != a21 && a12 != a22 {
                continue;
            }
            let atoms = if a12 == a21 {
                [a11, a12, a22]
            } else if a12 == a22 {
                [a11, a12, a21]
            } else if a11 == a21 {
                [a12, a11, a22]
            } else {
                [a12, a11, a21]
            };
            let oi = bond_orders.get(i).copied().unwrap_or(1);
            let oj = bond_orders.get(j).copied().unwrap_or(1);
            let triple = oi == 3
                || oj == 3
                || (oi == 2 && oj == 2 && degrees.get(atoms[1]).copied().unwrap_or(0) == 2);
            out.push(Angle { atoms, triple });
        }
    }
    out
}
