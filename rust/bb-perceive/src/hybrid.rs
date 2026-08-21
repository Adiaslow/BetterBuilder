//! Conjugation and hybridization.
//!
//! Hybridization reaches further than it looks: UFF atom typing keys on it, so it feeds bond
//! lengths and angles in `setTopolBounds`, and the UFF inversion terms apply only at SP2 centres.
//! It depends on conjugation, which depends in turn on `countAtomElec` — the same electron count
//! the aromaticity model uses.
//!
//! Computed over heavy atoms only. RDKit sanitizes before `AddHs`, so appended hydrogens never
//! have hybridization assigned and stay unspecified.

use crate::hydrogens::{contribution, Counts};
use crate::valence;
use crate::SmilesGraph;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hybridization {
    Unspecified = 0,
    S = 1,
    Sp = 2,
    Sp2 = 3,
    Sp3 = 4,
    Sp3d = 5,
    Sp3d2 = 6,
}

struct Ctx<'a> {
    g: &'a SmilesGraph,
    counts: &'a Counts,
    /// bond valence contributions, per bond
    contrib: Vec<f64>,
    /// neighbours as (atom, bond index)
    adj: Vec<Vec<(usize, usize)>>,
}

impl Ctx<'_> {
    fn total_degree(&self, i: usize) -> usize {
        self.adj[i].len() + self.counts.total_num_hs[i].max(0) as usize
    }

    /// The degree `numBondsPlusLonePairs` counts (ConjugHybrid.cpp): total degree less any bond this
    /// atom *donates* as a dative bond — RDKit decrements those (`isDative && idx != endAtomIdx`)
    /// before adding lone pairs, since the donated pair became the bond and is no longer a free pair.
    fn bonds_plus_lp_degree(&self, i: usize) -> i32 {
        let donated = self.adj[i]
            .iter()
            .filter(|&&(_, bi)| self.g.dative_donor[bi] == Some(i))
            .count() as i32;
        self.total_degree(i) as i32 - donated
    }

    fn total_valence(&self, i: usize) -> i32 {
        self.counts.explicit_valence[i] + self.counts.implicit_valence[i]
    }

    fn radicals(&self, i: usize) -> i32 {
        i32::from(self.counts.radicals[i])
    }

    /// Electrons available to the pi system — RDKit's `countAtomElec`.
    fn count_atom_elec(&self, i: usize) -> i32 {
        let z = self.g.atoms[i].atomic_number;
        let dv = valence::valence_list(z)
            .and_then(|v| v.first())
            .map_or(-1, |&v| i32::from(v));
        if dv <= 1 {
            return -1;
        }
        let degree = self.total_degree(i) as i32;
        if degree > 3 {
            return -1;
        }
        let nouter = valence::n_outer_elecs(z).map_or(0, i32::from);
        let mut nlp = (nouter - dv).max(0);
        nlp = (nlp - i32::from(self.g.atoms[i].charge)).max(0);
        // electrons available for donation into the pi system, RDKit `countAtomElec`
        // (Aromaticity.cpp): `(dv − degree) + nlp − nRadicals`. A radical electron is not available
        // to conjugate, so it is subtracted — without this a benzyl-type radical carbon reads as a
        // conjugation candidate and is wrongly hybridised SP2.
        let mut res = (dv - degree) + nlp - self.radicals(i);
        if res > 1 {
            // more than one unsaturation means a triple bond or higher: only one electron counts
            let unsat = self.counts.explicit_valence[i] - self.adj[i].len() as i32;
            if unsat > 1 {
                res = 1;
            }
        }
        res
    }

    fn is_conjug_cand(&self, i: usize) -> bool {
        let z = self.g.atoms[i].atomic_number;
        let vals = valence::valence_list(z).unwrap_or(&[]);
        let first = vals.first().copied().unwrap_or(-1);
        if self.g.atoms[i].charge == 0
            && first >= 0
            && self.total_valence(i) > i32::from(first)
        {
            return false;
        }
        let nouter = valence::n_outer_elecs(z).map_or(0, i32::from);
        ((z <= 10) || (nouter != 5 && nouter != 6) || (nouter == 6 && self.total_degree(i) < 2))
            && self.count_atom_elec(i) > 0
    }
}

/// Conjugated flag per bond, and hybridization per heavy atom.
pub fn perceive(g: &SmilesGraph, counts: &Counts) -> (Vec<bool>, Vec<Hybridization>) {
    let n = g.atoms.len();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    let mut contrib = Vec::with_capacity(g.bonds.len());
    for (bi, &(a, b)) in g.bonds.iter().enumerate() {
        let ar = g.atoms[a].aromatic_as_written && g.atoms[b].aromatic_as_written;
        contrib.push(contribution(g.bond_kinds[bi], ar));
        adj[a].push((b, bi));
        adj[b].push((a, bi));
    }
    let ctx = Ctx { g, counts, contrib, adj };

    // every aromatic bond is conjugated to begin with
    let mut conj: Vec<bool> = (0..g.bonds.len())
        .map(|bi| {
            let (a, b) = g.bonds[bi];
            g.atoms[a].aromatic_as_written && g.atoms[b].aromatic_as_written
        })
        .collect();

    for at in 0..n {
        if !ctx.is_conjug_cand(at) {
            continue;
        }
        let sbo = ctx.adj[at].len() + counts.total_num_hs[at].max(0) as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        for &(o1, b1) in &ctx.adj[at] {
            if ctx.contrib[b1] < 1.5 || !ctx.is_conjug_cand(o1) {
                continue;
            }
            for &(o2, b2) in &ctx.adj[at] {
                if b1 == b2 {
                    continue;
                }
                let sbo2 = ctx.adj[o2].len() + counts.total_num_hs[o2].max(0) as usize;
                if sbo2 > 3 {
                    continue;
                }
                if ctx.is_conjug_cand(o2) {
                    conj[b1] = true;
                    conj[b2] = true;
                }
            }
        }
    }

    let hyb = (0..n)
        .map(|i| {
            let z = g.atoms[i].atomic_number;
            if z == 0 {
                return Hybridization::Unspecified;
            }
            let deg = ctx.total_degree(i) as i32;
            // `numBondsPlusLonePairs` counts the dative-decremented degree; the norbs==4 SP2/SP3
            // decision below uses the raw total degree (RDKit's `getTotalDegree()`), so keep both.
            let lp_deg = ctx.bonds_plus_lp_degree(i);
            let norbs = if z < 89 {
                let nouter = valence::n_outer_elecs(z).map_or(0, i32::from);
                if z <= 1 {
                    lp_deg
                } else {
                    let tv = ctx.total_valence(i);
                    let chg = i32::from(g.atoms[i].charge);
                    let free = nouter - (tv + chg);
                    // RDKit `numBondsPlusLonePairs` (ConjugHybrid.cpp): norbs = degree + lone pairs.
                    // Below an octet it folds in radical electrons (radicals from `assignRadicals`);
                    // at or above an octet the radical term drops out.
                    if tv + nouter - chg < 8 {
                        let n_rad = ctx.radicals(i);
                        lp_deg + (free - n_rad) / 2 + n_rad
                    } else {
                        lp_deg + free / 2
                    }
                }
            } else {
                deg
            };
            match norbs {
                0 | 1 => Hybridization::S,
                2 => Hybridization::Sp,
                3 => Hybridization::Sp2,
                4 => {
                    let conjugated = ctx.adj[i].iter().any(|&(_, b)| conj[b]);
                    if deg > 3 || !conjugated {
                        Hybridization::Sp3
                    } else {
                        Hybridization::Sp2
                    }
                }
                5 => Hybridization::Sp3d,
                6 => Hybridization::Sp3d2,
                _ => Hybridization::Unspecified,
            }
        })
        .collect();

    (conj, hyb)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hyb_of(smiles: &str) -> Vec<Hybridization> {
        let g = crate::parse(smiles).unwrap();
        let c = crate::hydrogens::counts(&g);
        perceive(&g, &c).1
    }

    #[test]
    fn benzyl_radical_carbon_is_sp3() {
        // The benzyl radical carbon [CH2] has one radical electron. RDKit assigns it SP3, not SP2:
        // numBondsPlusLonePairs folds the radical in (norbs=4), and countAtomElec subtracts it so the
        // carbon is NOT a conjugation candidate — both radical terms are needed to reach SP3.
        let h = hyb_of("[CH2]c1ccccc1");
        assert_eq!(h[0], Hybridization::Sp3, "benzyl CH2 radical must be SP3");
        // the aromatic ring carbons stay SP2
        assert_eq!(h[1], Hybridization::Sp2);
    }

    #[test]
    fn closed_shell_hybridization_unaffected() {
        // sanity anchors that the radical terms are no-ops without radicals
        assert_eq!(hyb_of("CC")[0], Hybridization::Sp3); // ethane C
        assert_eq!(hyb_of("C=C")[0], Hybridization::Sp2); // ethylene C
        assert_eq!(hyb_of("C#C")[0], Hybridization::Sp); // acetylene C
        assert_eq!(hyb_of("CC(=O)O")[2], Hybridization::Sp2); // carboxyl C=O and both O sp2-ish
    }
}
