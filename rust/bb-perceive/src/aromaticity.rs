//! Aromaticity perception — a faithful port of RDKit's default model (`GraphMol/Aromaticity.cpp`,
//! `aromaticityHelper(mol, srings, 0, 0, includeFused=true)`). RDKit re-perceives aromaticity from the
//! kekulized structure regardless of how the SMILES wrote it (lowercase vs Kekulé `C=C`), so this
//! reproduces that: given the kekulized bond orders, ring set, and per-atom perception data, it returns
//! which atoms and bonds are aromatic.
//!
//! The electron-donor typing keys on `countAtomElec` (the same pi-electron count as [`crate::hybrid`]);
//! the Hückel rule is applied over every fused-ring subsystem combinatorially, exactly as RDKit does.

use crate::valence;

/// Per-atom inputs the perception needs — all already computed elsewhere in perception, in atom order.
pub struct AromInput<'a> {
    pub atomic_number: &'a [u8],
    pub charge: &'a [i8],
    pub radicals: &'a [u8],
    pub total_num_hs: &'a [i32],
    /// Kekulized bond orders (1 single, 2 double, 3 triple) per bond, and the bond endpoints.
    pub bonds: &'a [(usize, usize)],
    pub bond_orders: &'a [u8],
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Edon {
    Vacant,
    One,
    Two,
    Any,
    None,
}

struct Mol<'a> {
    n: usize,
    inp: &'a AromInput<'a>,
    adj: Vec<Vec<usize>>,       // atom -> bond indices
    bond_in_ring: Vec<bool>,    // per bond
    degree: Vec<usize>,         // heavy-atom degree (bond count)
}

impl<'a> Mol<'a> {
    fn new(inp: &'a AromInput<'a>, n: usize, rings: &[Vec<usize>]) -> Self {
        let mut adj = vec![Vec::new(); n];
        for (bi, &(a, b)) in inp.bonds.iter().enumerate() {
            adj[a].push(bi);
            adj[b].push(bi);
        }
        // a bond is "in a ring" iff both endpoints appear consecutively in some SSSR ring
        let mut bond_in_ring = vec![false; inp.bonds.len()];
        let ring_bond: std::collections::HashSet<(usize, usize)> = rings
            .iter()
            .flat_map(|r| r.iter().enumerate().map(move |(i, &a)| {
                let b = r[(i + 1) % r.len()];
                if a < b { (a, b) } else { (b, a) }
            }))
            .collect();
        for (bi, &(a, b)) in inp.bonds.iter().enumerate() {
            let key = if a < b { (a, b) } else { (b, a) };
            bond_in_ring[bi] = ring_bond.contains(&key);
        }
        let degree = adj.iter().map(Vec::len).collect();
        Mol { n, inp, adj, bond_in_ring, degree }
    }

    fn other(&self, bi: usize, at: usize) -> usize {
        let (a, b) = self.inp.bonds[bi];
        if a == at { b } else { a }
    }

    fn total_degree(&self, i: usize) -> i32 {
        self.degree[i] as i32 + self.inp.total_num_hs[i].max(0)
    }

    fn total_valence(&self, i: usize) -> i32 {
        // sum of bond orders + hydrogens (kekulized: integer orders)
        let bo: i32 = self.adj[i].iter().map(|&bi| i32::from(self.inp.bond_orders[bi])).sum();
        bo + self.inp.total_num_hs[i].max(0)
    }

    /// RDKit `countAtomElec` (Aromaticity.cpp) — pi electrons available for donation.
    fn count_atom_elec(&self, i: usize) -> i32 {
        let z = self.inp.atomic_number[i];
        let dv = valence::valence_list(z).and_then(|v| v.first()).map_or(-1, |&v| i32::from(v));
        if dv <= 1 {
            return -1;
        }
        let degree = self.total_degree(i);
        if degree > 3 {
            return -1;
        }
        let nouter = valence::n_outer_elecs(z).map_or(0, i32::from);
        let mut nlp = (nouter - dv).max(0);
        nlp = (nlp - i32::from(self.inp.charge[i])).max(0);
        let mut res = (dv - degree) + nlp - i32::from(self.inp.radicals[i]);
        if res > 1 {
            // more than one unsaturation → triple bond, only one electron counts
            let explicit: i32 = self.adj[i].iter().map(|&bi| i32::from(self.inp.bond_orders[bi])).sum::<i32>()
                + self.inp.total_num_hs[i].max(0);
            let n_unsat = explicit - self.degree[i] as i32;
            if n_unsat > 1 {
                res = 1;
            }
        }
        res
    }

    fn more_electro_negative(&self, z1: u8, z2: u8) -> bool {
        let ne1 = valence::n_outer_elecs(z1).unwrap_or(0);
        let ne2 = valence::n_outer_elecs(z2).unwrap_or(0);
        ne1 > ne2 || (ne1 == ne2 && z1 < z2)
    }

    fn incident_noncyclic_multiple(&self, at: usize) -> Option<usize> {
        for &bi in &self.adj[at] {
            if !self.bond_in_ring[bi] && self.inp.bond_orders[bi] >= 2 {
                return Some(self.other(bi, at));
            }
        }
        None
    }

    fn incident_cyclic_multiple(&self, at: usize) -> bool {
        self.adj[at].iter().any(|&bi| self.bond_in_ring[bi] && self.inp.bond_orders[bi] >= 2)
    }

    fn incident_multiple(&self, at: usize) -> bool {
        // RDKit: explicitValence != (degree adjusted for zero-order bonds). With integer kekulized
        // orders, this is simply "has a bond of order >= 2".
        self.adj[at].iter().any(|&bi| self.inp.bond_orders[bi] >= 2)
    }

    /// RDKit `getAtomDonorTypeArom` (exocyclicBondsStealElectrons = true).
    fn donor_type(&self, at: usize) -> Edon {
        let mut nelec = self.count_atom_elec(at);
        if nelec < 0 {
            return Edon::None;
        }
        if nelec == 0 {
            return if self.incident_noncyclic_multiple(at).is_some() {
                Edon::Vacant
            } else if self.incident_cyclic_multiple(at) {
                Edon::One
            } else {
                Edon::None
            };
        }
        if nelec == 1 {
            if let Some(who) = self.incident_noncyclic_multiple(at) {
                return if self.more_electro_negative(self.inp.atomic_number[who], self.inp.atomic_number[at]) {
                    Edon::Vacant
                } else {
                    Edon::One
                };
            }
            if self.incident_multiple(at) {
                return Edon::One;
            }
            if self.inp.charge[at] == 1 {
                return Edon::Vacant;
            }
            return Edon::None;
        }
        // nelec >= 2
        if let Some(who) = self.incident_noncyclic_multiple(at) {
            if self.more_electro_negative(self.inp.atomic_number[who], self.inp.atomic_number[at]) {
                nelec -= 1;
            }
        }
        if nelec % 2 == 1 { Edon::One } else { Edon::Two }
    }

    /// RDKit `isAtomCandForArom` (default flags: allowThirdRow, allowTripleBonds, allowHigherExceptions).
    fn is_cand(&self, at: usize, edon: Edon) -> bool {
        let z = self.inp.atomic_number[at];
        if z > 18 && z != 34 && z != 52 {
            return false;
        }
        if matches!(edon, Edon::None) {
            return false;
        }
        // not in default valence state → shut out
        let def_val = valence::valence_list(z).and_then(|v| v.first()).map_or(-1, |&v| i32::from(v));
        if def_val > 0 {
            let eff = (i32::from(z) - i32::from(self.inp.charge[at])).clamp(0, 118) as u8;
            let eff_dv = valence::valence_list(eff).and_then(|v| v.first()).map_or(-1, |&v| i32::from(v));
            if self.total_valence(at) > eff_dv {
                return false;
            }
        }
        // radical heteroatoms / charged radical carbons disqualify
        if self.inp.radicals[at] != 0 && (z != 6 || self.inp.charge[at] != 0) {
            return false;
        }
        // disallow more than one multiple bond (unless one is a triple; allowTripleBonds=true)
        let explicit: i32 = self.adj[at].iter().map(|&bi| i32::from(self.inp.bond_orders[bi])).sum::<i32>()
            + self.inp.total_num_hs[at].max(0);
        let n_unsat = explicit - self.degree[at] as i32;
        if n_unsat > 1 {
            let n_mult = self.adj[at].iter().filter(|&&bi| self.inp.bond_orders[bi] >= 2).count();
            // allowTripleBonds: a single triple bond (n_unsat 2, n_mult 1) is allowed
            if !(n_mult == 1 && n_unsat == 2) {
                return false;
            }
        }
        true
    }
}

fn min_max_elecs(d: Edon) -> (i32, i32) {
    match d {
        Edon::Any => (1, 2),
        Edon::One => (1, 1),
        Edon::Two => (2, 2),
        Edon::None | Edon::Vacant => (0, 0),
    }
}

/// RDKit `applyHuckel`: does this atom set satisfy 4n+2 given each atom's donor range?
fn apply_huckel(atoms: &[usize], edon: &[Edon]) -> bool {
    let (mut rlw, mut rup) = (0i32, 0i32);
    let mut n_any = 0;
    for &idx in atoms {
        if edon[idx] == Edon::Any {
            n_any += 1;
            if n_any > 1 {
                return false;
            }
        }
        let (lw, up) = min_max_elecs(edon[idx]);
        rlw += lw;
        rup += up;
    }
    if rup >= 6 {
        (rlw..=rup).any(|rie| (rie - 2) % 4 == 0)
    } else {
        rup == 2
    }
}

/// RDKit `nextCombination`.
fn next_combination(comb: &mut [i32], tot: i32) -> i32 {
    let nelem = comb.len() as i32;
    let mut celem = nelem - 1;
    while comb[celem as usize] == tot - nelem + celem {
        celem -= 1;
        if celem < 0 {
            return -1;
        }
    }
    comb[celem as usize] += 1;
    for i in (celem as usize + 1)..comb.len() {
        comb[i] = comb[i - 1] + 1;
    }
    celem
}

/// Perceive aromaticity. Returns `(atom_aromatic, bond_aromatic)`.
pub fn perceive(inp: &AromInput, rings: &[Vec<usize>]) -> (Vec<bool>, Vec<bool>) {
    let n = inp.atomic_number.len();
    let mol = Mol::new(inp, n, rings);

    // donor type + candidacy per atom, and the candidate rings (all atoms aromatic-candidates)
    let mut edon = vec![Edon::None; n];
    let mut cand = vec![false; n];
    let mut seen = vec![false; n];
    let mut c_rings: Vec<Vec<usize>> = Vec::new();
    for ring in rings {
        let mut all_arom = true;
        for &idx in ring {
            if !seen[idx] {
                seen[idx] = true;
                edon[idx] = mol.donor_type(idx);
                cand[idx] = mol.is_cand(idx, edon[idx]);
            }
            if !cand[idx] {
                all_arom = false;
            }
        }
        if all_arom {
            c_rings.push(ring.clone());
        }
    }

    // rings as bond-id sets
    let bond_id = |a: usize, b: usize| -> usize {
        mol.adj[a].iter().copied().find(|&bi| mol.other(bi, a) == b).expect("ring bond exists")
    };
    let brings: Vec<Vec<usize>> = c_rings
        .iter()
        .map(|r| (0..r.len()).map(|i| bond_id(r[i], r[(i + 1) % r.len()])).collect())
        .collect();

    // Ring neighbor map — RDKit RingUtils::makeRingNeighborMap(brings, neighMap, 24, 1): two rings are
    // neighbors only if they share EXACTLY 1 bond (maxOverlapSize = 1) and neither has more than
    // `MAX_FUSED_RING_SIZE` bonds. This is what stops, e.g., a porphyrin's 16-membered macrocycle
    // (which shares 2 bonds with each pyrrole) from making the whole system one fused set — RDKit then
    // evaluates only the single rings, so a non-aromatic pyrrole stays non-aromatic.
    const MAX_FUSED_RING_SIZE: usize = 24;
    let ncr = c_rings.len();
    let mut neigh: Vec<Vec<usize>> = vec![Vec::new(); ncr];
    for i in 0..ncr {
        if brings[i].len() > MAX_FUSED_RING_SIZE {
            continue;
        }
        let si: std::collections::HashSet<usize> = brings[i].iter().copied().collect();
        for j in (i + 1)..ncr {
            if brings[j].len() > MAX_FUSED_RING_SIZE {
                continue;
            }
            let shared = brings[j].iter().filter(|b| si.contains(b)).count();
            if shared > 0 && shared <= 1 {
                neigh[i].push(j);
                neigh[j].push(i);
            }
        }
    }

    let mut atom_arom = vec![false; n];
    let mut bond_arom = vec![false; inp.bonds.len()];

    // pick each fused system, apply Hückel to its subsystems
    let mut fus_done = vec![false; ncr];
    let mut curr = 0usize;
    while curr < ncr {
        let mut fused = Vec::new();
        pick_fused(curr, &neigh, &mut fused, &mut fus_done);
        apply_huckel_to_fused(&mol, &c_rings, &brings, &fused, &edon, &neigh, &mut atom_arom, &mut bond_arom);
        match (0..ncr).find(|&r| !fus_done[r]) {
            Some(r) => curr = r,
            None => break,
        }
    }

    (atom_arom, bond_arom)
}

fn pick_fused(curr: usize, neigh: &[Vec<usize>], res: &mut Vec<usize>, done: &mut [bool]) {
    done[curr] = true;
    res.push(curr);
    for &nb in &neigh[curr] {
        if !done[nb] {
            pick_fused(nb, neigh, res, done);
        }
    }
}

/// RingUtils::checkFused — are exactly these ring ids a single fused system?
fn check_fused(rids: &[usize], neigh: &[Vec<usize>]) -> bool {
    let nrings = neigh.len();
    let mut done = vec![false; nrings];
    for (r, d) in done.iter_mut().enumerate() {
        // mark every ring outside `rids` as done, so pick_fused only walks the candidate set
        *d = !rids.contains(&r);
    }
    let mut fused = Vec::new();
    pick_fused(rids[0], neigh, &mut fused, &mut done);
    fused.len() == rids.len()
}

#[allow(clippy::too_many_arguments)]
fn apply_huckel_to_fused(
    mol: &Mol,
    srings: &[Vec<usize>],
    brings: &[Vec<usize>],
    fused: &[usize],
    edon: &[Edon],
    neigh: &[Vec<usize>],
    atom_arom: &mut [bool],
    bond_arom: &mut [bool],
) {
    let nrings = fused.len() as i32;
    let n_ring_bonds: usize = {
        let mut set = std::collections::HashSet::new();
        for &ridx in fused {
            set.extend(brings[ridx].iter().copied());
        }
        set.len()
    };
    let mut done_bonds: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let max_fused = 6i32;

    let mut cur_size = 0i32;
    let mut comb: Vec<i32> = Vec::new();
    let mut pos = -1i32;
    loop {
        if pos == -1 {
            if cur_size == 2 && nrings > 300 {
                break;
            }
            cur_size += 1;
            if cur_size > nrings.min(max_fused) || done_bonds.len() >= n_ring_bonds {
                break;
            }
            comb = (0..cur_size).collect();
            pos = 0;
        } else {
            pos = next_combination(&mut comb, nrings);
        }
        if pos == -1 {
            continue;
        }
        let cur_rs: Vec<usize> = comb.iter().map(|&i| fused[i as usize]).collect();
        if !check_fused(&cur_rs, neigh) {
            continue;
        }
        // atoms present in one or two of the subsystem rings
        let mut counts = vec![0i32; mol.n];
        for &ridx in &cur_rs {
            for &a in &srings[ridx] {
                counts[a] += 1;
            }
        }
        let unon: Vec<usize> = (0..mol.n).filter(|&i| counts[i] == 1 || counts[i] == 2).collect();
        if apply_huckel(&unon, edon) {
            mark_arom(brings, &cur_rs, &mut done_bonds, mol, atom_arom, bond_arom);
        }
    }
}

/// RDKit `markAtomsBondsArom`: only bonds appearing once across the fused subsystem (the perimeter)
/// are marked aromatic — and their atoms with them.
fn mark_arom(
    brings: &[Vec<usize>],
    ring_ids: &[usize],
    done_bonds: &mut std::collections::HashSet<usize>,
    mol: &Mol,
    atom_arom: &mut [bool],
    bond_arom: &mut [bool],
) {
    let mut bnd_cnt: std::collections::HashMap<usize, i32> = std::collections::HashMap::new();
    for &ri in ring_ids {
        for &bi in &brings[ri] {
            *bnd_cnt.entry(bi).or_insert(0) += 1;
        }
    }
    for (&bi, &cnt) in &bnd_cnt {
        if cnt == 1 {
            bond_arom[bi] = true;
            let (a, b) = mol.inp.bonds[bi];
            atom_arom[a] = true;
            atom_arom[b] = true;
            done_bonds.insert(bi);
        }
    }
}
