//! Chirality and stereo constraints for the distance-geometry problem.
//!
//! Ports RDKit's `EmbeddingOps::findDoubleBonds` and `findChiralSets` (`Embedder.cpp`), which
//! produce the `MoleculeSpec` fields the acceptance checks read: double-bond substituent triples,
//! cis/trans stereo bonds, chiral-volume sets, and untagged tetrahedral centers.

use bb_perceive::chirality::ChiralTag;
use bb_perceive::smarts_match::{Perceived, Stereo};

/// A chiral-volume set: `center`, four substituent atoms, signed volume bounds, and the fused-small-
/// ring flag. Matches RDKit's `DistGeom::ChiralSet` layout (`d_idx0`, `d_idx1..4`, vol bounds).
#[derive(Debug, Clone, PartialEq)]
pub struct ChiralSet {
    pub center: u32,
    pub atoms: [u32; 4],
    pub vol_lo: f32,
    pub vol_hi: f32,
    pub fused_small_rings: bool,
}

/// Reproduce `findChiralSets`: tagged stereocenters become signed-volume energy terms
/// (`chiral_sets`), untagged C/N degree-4 centers in ≥2 rings become non-degeneracy checks
/// (`tetrahedral_centers`, zero volume). Neighbours are taken in bond order (RDKit `getAtomBonds`),
/// which is gated to match RDKit; a 3-coordinate center appends itself as the fourth point.
pub fn find_chiral_sets(mol: &Perceived) -> (Vec<ChiralSet>, Vec<ChiralSet>) {
    let n = mol.atomic_numbers.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b) in &mol.bonds {
        adj[a].push(b);
        adj[b].push(a);
    }
    let num_atom_rings = |atom: usize| mol.rings.iter().filter(|r| r.contains(&atom)).count();
    let in_ring_of_size =
        |atom: usize, sz: usize| mol.rings.iter().any(|r| r.len() == sz && r.contains(&atom));

    let mut chiral = Vec::new();
    let mut tetra = Vec::new();
    for (x, neighbours) in adj.iter().enumerate() {
        if mol.atomic_numbers[x] == 1 {
            continue; // skip hydrogens
        }
        let tag = mol.chiral_tags[x];
        let z = mol.atomic_numbers[x];
        let degree = neighbours.len();
        let tagged = matches!(tag, ChiralTag::Cw | ChiralTag::Ccw);
        if !(tagged || ((z == 6 || z == 7) && degree == 4)) {
            continue;
        }
        let mut nbrs = neighbours.clone();
        let vol_lower = if nbrs.len() < 4 {
            nbrs.push(x); // include the center as the fourth point
            2.0f32
        } else {
            5.0f32
        };
        let vol_upper = 100.0f32;
        let fused = mol
            .rings
            .iter()
            .filter(|r| r.len() < 5 && r.contains(&x))
            .count()
            > 1;
        let atoms = [
            nbrs[0] as u32,
            nbrs[1] as u32,
            nbrs[2] as u32,
            nbrs[3] as u32,
        ];
        match tag {
            ChiralTag::Ccw => chiral.push(ChiralSet {
                center: x as u32,
                atoms,
                vol_lo: vol_lower,
                vol_hi: vol_upper,
                fused_small_rings: fused,
            }),
            ChiralTag::Cw => chiral.push(ChiralSet {
                center: x as u32,
                atoms,
                vol_lo: -vol_upper,
                vol_hi: -vol_lower,
                fused_small_rings: fused,
            }),
            ChiralTag::None => {
                // untagged C/N degree-4: a non-degeneracy check, but only for atoms in >=2 rings
                // and not in a 3-ring
                if !(num_atom_rings(x) < 2 || in_ring_of_size(x, 3)) {
                    tetra.push(ChiralSet {
                        center: x as u32,
                        atoms,
                        vol_lo: 0.0,
                        vol_hi: 0.0,
                        fused_small_rings: fused,
                    });
                }
            }
        }
    }
    (chiral, tetra)
}

/// A `(nbr, dbAtom, otherDbAtom)` triple for the double-bond linear-arrangement check.
pub type DoubleBondEnd = [u32; 3];

/// A stereo double bond: controlling atoms `[stereoAtom0, begin, end, stereoAtom1]` and `sign`
/// (+1 trans/E, −1 cis/Z).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoDoubleBond {
    pub atoms: [u32; 4],
    pub sign: i8,
}

fn adjacency(mol: &Perceived) -> Vec<Vec<(usize, usize)>> {
    let n = mol.atomic_numbers.len();
    let mut adj = vec![Vec::new(); n];
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        adj[a].push((bi, b));
        adj[b].push((bi, a));
    }
    adj
}

/// Reproduce `findDoubleBonds`: the substituent triples and stereo double bonds.
///
/// `doubleBondEnds` is pure topology. `stereoDoubleBonds` uses the perceived double-bond stereo
/// ([`bb_perceive::stereo`]); its `[stereoAtom0, begin, end, stereoAtom1]`/`sign` encoding picks the
/// directional-bond neighbour as each end's reference atom, whereas RDKit's legacy
/// `assignBondStereoCodes` picks the higher-CIP-rank neighbour. The two choices coincide unless a
/// double-bond carbon carries two heavy substituents and the directional bond is on the lower-CIP
/// one (e.g. `/C(F)=C/`, where F outranks the directional carbon): there the reference and the Z/E
/// label both flip. The two in-plane substituents lie ~180° apart across the C=C axis, so switching
/// the reference shifts the enforced dihedral by π and the compensating sign flip leaves the cis/trans
/// acceptance check (`bb_embed::checks::double_bond_stereo_ok`) accepting exactly the same conformers
/// — and the bounds matrix is identical either way (the same representation-independence as set14,
/// gated by `validation/parity/gate_spec.py`).
pub fn find_double_bonds(mol: &Perceived) -> (Vec<DoubleBondEnd>, Vec<StereoDoubleBond>) {
    let adj = adjacency(mol);
    let mut ends = Vec::new();
    let mut stereo = Vec::new();
    for (bi, &(begin, end)) in mol.bonds.iter().enumerate() {
        if mol.bond_order[bi] != 2 {
            continue;
        }
        for &atm in &[begin, end] {
            let degree = adj[atm].len();
            if degree < 2 {
                continue;
            }
            let oatm = if atm == begin { end } else { begin };
            for &(obnd, nbr) in &adj[atm] {
                if nbr == oatm {
                    continue;
                }
                // a degree-2 double-bond atom only counts its single-bond neighbours
                if mol.bond_order[obnd] != 1 && degree == 2 {
                    continue;
                }
                ends.push([nbr as u32, atm as u32, oatm as u32]);
            }
        }
        // stereo: definite (> Any) double bonds
        if matches!(
            mol.bond_stereo[bi],
            Stereo::Z | Stereo::E | Stereo::Cis | Stereo::Trans
        ) {
            let sa = mol.stereo_atoms[bi];
            let sign = if matches!(mol.bond_stereo[bi], Stereo::Cis | Stereo::Z) {
                -1
            } else {
                1
            };
            stereo.push(StereoDoubleBond {
                atoms: [sa[0] as u32, begin as u32, end as u32, sa[1] as u32],
                sign,
            });
        }
    }
    (ends, stereo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bb_perceive::smarts_match::perceive;

    #[test]
    fn double_bond_ends_of_a_simple_alkene() {
        // C/C=C/C  (2-butene): the central C=C, each end a CH3 (degree... post-AddHs the end C is
        // degree 4). Every non-partner neighbour of each double-bond atom is a substituent triple.
        let mol = perceive("CC=CC").expect("perceives");
        let (ends, _) = find_double_bonds(&mol);
        // each double-bond carbon has the methyl C plus (post-AddHs) one H as substituents
        assert!(ends.iter().any(|e| e[1] == 1 && e[2] == 2)); // begin atom 1, partner 2
        assert!(ends.iter().any(|e| e[1] == 2 && e[2] == 1)); // end atom 2, partner 1
    }
}
