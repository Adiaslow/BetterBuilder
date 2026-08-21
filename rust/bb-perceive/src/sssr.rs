//! Ring perception reproducing what `SANITIZE_SYMMRINGS` leaves on a molecule.
//!
//! The result is a minimum cycle basis together with the symmetry-equivalent rings RDKit restores,
//! because that set is what every downstream query reads: ring-size SMARTS such as `r{9-}`, the
//! ring-atom union behind `count_exo_rotatable`, and largest-ring selection for the pinned core.
//!
//! Which representative survives among equal-size alternates is decided by traversal, not by a
//! closed-form rule, so the mechanism is reproduced rather than a criterion invented for it:
//! degree-2 chains contribute one search root each, rings are kept the first time their atom set is
//! seen, nodes that discover the same ring are re-searched with each other's bonds removed, and the
//! surviving set is pruned by bond coverage. See `rings.rs` for an independent minimum-cycle-basis
//! implementation used as a cross-check.

use std::collections::{HashSet, VecDeque};

const WHITE: u8 = 0;
const GRAY: u8 = 1;
const BLACK: u8 = 2;

pub struct Rings {
    n: usize,
    /// `adj[a]` is `(neighbour, bond index)`.
    adj: Vec<Vec<(usize, usize)>>,
    bonds: Vec<(usize, usize)>,
}

/// Put a ring walk in RDKit's canonical rotation, in place.
///
/// The smallest atom index leads, and the walk then runs towards whichever of its two ring
/// neighbours is smaller. `FindRings.cpp`'s `normalize_ring` applies this to every ring before the
/// extra-ring cleanup, so it is the form rings are both compared and stored in.
fn normalize_ring(ring: &mut [usize]) {
    if ring.len() < 3 {
        return;
    }
    let start = ring
        .iter()
        .enumerate()
        .min_by_key(|&(_, v)| *v)
        .map(|(i, _)| i)
        .expect("ring is non-empty");
    ring.rotate_left(start);
    if ring[ring.len() - 1] < ring[1] {
        // reverse the walk while leaving the smallest atom in front
        ring[1..].reverse();
    }
}

/// A ring's identity is its atom set, so the same cycle found by different routes dedups.
fn invariant(ring: &[usize]) -> Vec<usize> {
    let mut v = ring.to_vec();
    v.sort_unstable();
    v
}

impl Rings {
    pub fn new(n_atoms: usize, bonds: &[(usize, usize)]) -> Self {
        let mut adj = vec![Vec::new(); n_atoms];
        for (bi, &(a, b)) in bonds.iter().enumerate() {
            adj[a].push((b, bi));
            adj[b].push((a, bi));
        }
        Self {
            n: n_atoms,
            adj,
            bonds: bonds.to_vec(),
        }
    }

    /// Deactivate every live bond at `cand`, updating degrees and queueing neighbours that drop
    /// to degree 1 or below.
    fn trim_bonds(
        &self,
        cand: usize,
        changed: &mut VecDeque<usize>,
        degrees: &mut [i32],
        active: &mut [bool],
    ) {
        for &(other, bi) in &self.adj[cand] {
            if !active[bi] {
                continue;
            }
            if degrees[other] <= 2 {
                changed.push_back(other);
            }
            active[bi] = false;
            degrees[other] -= 1;
            degrees[cand] -= 1;
        }
    }

    /// Mark degree-2 atoms reachable from `root` through degree-2 atoms: one root per chain is
    /// enough to find the same rings.
    fn mark_useless_d2s(&self, root: usize, forb: &mut [bool], degrees: &[i32], active: &[bool]) {
        let mut stack = vec![root];
        while let Some(a) = stack.pop() {
            for &(other, bi) in &self.adj[a] {
                if active[bi] && !forb[other] && degrees[other] == 2 {
                    forb[other] = true;
                    stack.push(other);
                }
            }
        }
    }

    fn pick_d2_nodes(&self, frag: &[usize], degrees: &[i32], active: &[bool]) -> Vec<usize> {
        let mut out = Vec::new();
        let mut forb = vec![false; self.n];
        loop {
            let root = frag
                .iter()
                .copied()
                .find(|&a| degrees[a] == 2 && !forb[a]);
            match root {
                None => break,
                Some(r) => {
                    out.push(r);
                    forb[r] = true;
                    self.mark_useless_d2s(r, &mut forb, degrees, active);
                }
            }
        }
        out
    }

    /// All rings of the smallest size through `root`, by BFS. Nodes already dequeued are closed to
    /// further ring closure, and the search stops once a larger ring appears.
    fn smallest_rings_bfs(
        &self,
        root: usize,
        active: &[bool],
        forbidden: Option<&[usize]>,
    ) -> Vec<Vec<usize>> {
        let mut done = vec![WHITE; self.n];
        if let Some(f) = forbidden {
            for &i in f {
                done[i] = BLACK;
            }
        }
        let mut parents = vec![usize::MAX; self.n];
        let mut depths = vec![0usize; self.n];
        let mut rings: Vec<Vec<usize>> = Vec::new();
        let mut q = VecDeque::from(vec![root]);
        let mut cur_size = usize::MAX;

        while let Some(curr) = q.pop_front() {
            done[curr] = BLACK;
            let depth = depths[curr] + 1;
            if depth > cur_size {
                break;
            }
            for &(nbr, bi) in &self.adj[curr] {
                if !active[bi] || done[nbr] == BLACK || parents[curr] == nbr {
                    continue;
                }
                if done[nbr] == WHITE {
                    parents[nbr] = curr;
                    done[nbr] = GRAY;
                    depths[nbr] = depth;
                    q.push_back(nbr);
                } else {
                    // reached by another path: stitch the two together
                    let mut ring = vec![nbr];
                    let mut p = parents[nbr];
                    while p != usize::MAX && p != root {
                        ring.push(p);
                        p = parents[p];
                    }
                    ring.insert(0, curr);
                    let mut p = parents[curr];
                    let mut ok = true;
                    while p != usize::MAX {
                        if ring.contains(&p) {
                            ok = false;
                            break;
                        }
                        ring.insert(0, p);
                        p = parents[p];
                    }
                    if !ok || ring.len() <= 1 {
                        continue;
                    }
                    if ring.len() <= cur_size {
                        cur_size = ring.len();
                        rings.push(ring);
                    } else {
                        return rings;
                    }
                }
            }
        }
        rings
    }

    /// Rings through the current degree-2 roots, plus a second pass for roots that found the same
    /// ring: each is re-searched with the others' bonds removed, which is what separates
    /// equal-size alternates.
    #[allow(clippy::too_many_arguments)]
    fn find_rings_d2_nodes(
        &self,
        res: &mut Vec<Vec<usize>>,
        invars: &mut HashSet<Vec<usize>>,
        d2nodes: &[usize],
        degrees: &mut [i32],
        active: &mut [bool],
    ) {
        let mut dup_cands: Vec<(Vec<usize>, Vec<usize>)> = Vec::new(); // (invariant, roots)
        let mut dup_map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();

        for &cand in d2nodes {
            let srings = self.smallest_rings_bfs(cand, active, None);
            for ring in &srings {
                let inv = invariant(ring);
                let slot = match dup_cands.iter().position(|(i, _)| *i == inv) {
                    Some(p) => p,
                    None => {
                        dup_cands.push((inv.clone(), Vec::new()));
                        dup_cands.len() - 1
                    }
                };
                if !invars.contains(&inv) {
                    res.push(ring.clone());
                    invars.insert(inv);
                } else {
                    for &other in &dup_cands[slot].1 {
                        dup_map.entry(cand).or_default().push(other);
                        dup_map.entry(other).or_default().push(cand);
                    }
                }
                dup_cands[slot].1.push(cand);
            }
            if srings.is_empty() {
                let mut changed = VecDeque::from(vec![cand]);
                while let Some(c) = changed.pop_front() {
                    self.trim_bonds(c, &mut changed, degrees, active);
                }
            }
        }

        // roots that found the same ring: re-search each with the others removed
        for (_, roots) in dup_cands.iter().filter(|(_, r)| r.len() > 1) {
            let mut nrings: Vec<Vec<usize>> = Vec::new();
            let mut min_size = usize::MAX;
            for &root in roots {
                let mut deg_copy = degrees.to_vec();
                let mut act_copy = active.to_vec();
                let mut changed = VecDeque::new();
                if let Some(others) = dup_map.get(&root) {
                    for &o in others {
                        self.trim_bonds(o, &mut changed, &mut deg_copy, &mut act_copy);
                    }
                }
                for ring in self.smallest_rings_bfs(root, &act_copy, None) {
                    min_size = min_size.min(ring.len());
                    nrings.push(ring);
                }
            }
            for ring in nrings.into_iter().filter(|r| r.len() == min_size) {
                let inv = invariant(&ring);
                if !invars.contains(&inv) {
                    res.push(ring);
                    invars.insert(inv);
                }
            }
        }
    }

    fn find_rings_d3_node(
        &self,
        res: &mut Vec<Vec<usize>>,
        invars: &mut HashSet<Vec<usize>>,
        cand: usize,
        active: &[bool],
    ) {
        for ring in self.smallest_rings_bfs(cand, active, None) {
            let inv = invariant(&ring);
            if !invars.contains(&inv) {
                res.push(ring);
                invars.insert(inv);
            }
        }
    }

    /// Prune to the rings that actually contribute bonds, keeping the discards as substitution
    /// candidates. Among equal-size rings the one overlapping most with what is already covered is
    /// considered first, which is what decides the survivor.
    fn remove_extra_rings(&self, mut cands: Vec<Vec<usize>>) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
        // `removeExtraRings` sorts the ring list itself, so both the rings it keeps and the extras
        // it sets aside come out size-ascending rather than in the order they were found.
        cands.sort_by_key(Vec::len);
        let order: Vec<usize> = (0..cands.len()).collect();

        let bond_sets: Vec<HashSet<usize>> =
            cands.iter().map(|r| self.ring_bond_set(r)).collect();

        let mut munion: HashSet<usize> = HashSet::new();
        let mut keep = vec![false; cands.len()];
        let mut avail = vec![true; cands.len()];

        for pos in 0..order.len() {
            let i = order[pos];
            if bond_sets[i].is_subset(&munion) {
                avail[i] = false;
            }
            if !avail[i] {
                continue;
            }
            munion.extend(&bond_sets[i]);
            keep[i] = true;

            let mut consider: Vec<usize> = order[pos + 1..]
                .iter()
                .copied()
                .take_while(|&j| cands[j].len() == cands[i].len())
                .filter(|&j| avail[j])
                .collect();
            while !consider.is_empty() {
                let best = consider
                    .iter()
                    .copied()
                    .max_by_key(|&j| bond_sets[j].intersection(&munion).count())
                    .expect("non-empty");
                consider.retain(|&j| j != best);
                if bond_sets[best].is_subset(&munion) {
                    avail[best] = false;
                } else {
                    keep[best] = true;
                    avail[best] = false;
                    munion.extend(&bond_sets[best]);
                }
            }
        }

        let mut kept = Vec::new();
        let mut extras = Vec::new();
        for (i, ring) in cands.into_iter().enumerate() {
            if keep[i] {
                kept.push(ring);
            } else {
                extras.push(ring);
            }
        }
        (kept, extras)
    }


    /// Bonds with both ends inside `frag`.
    fn fragment_bond_count(&self, frag: &[usize]) -> usize {
        let members: HashSet<usize> = frag.iter().copied().collect();
        self.bonds
            .iter()
            .filter(|(a, b)| members.contains(a) && members.contains(b))
            .count()
    }

    fn ring_bond_set(&self, ring: &[usize]) -> HashSet<usize> {
        let atoms: HashSet<usize> = ring.iter().copied().collect();
        let mut out = HashSet::new();
        for (bi, &(a, b)) in self.bonds.iter().enumerate() {
            if atoms.contains(&a) && atoms.contains(&b) {
                out.insert(bi);
            }
        }
        // only bonds that are part of the cycle itself
        if out.len() > ring.len() {
            out.retain(|&bi| {
                let (a, b) = self.bonds[bi];
                ring_adjacent(ring, a, b)
            });
        }
        out
    }

    fn fragments(&self) -> Vec<Vec<usize>> {
        let mut seen = vec![false; self.n];
        let mut out = Vec::new();
        for start in 0..self.n {
            if seen[start] {
                continue;
            }
            let mut frag = Vec::new();
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(a) = stack.pop() {
                frag.push(a);
                for &(nb, _) in &self.adj[a] {
                    if !seen[nb] {
                        seen[nb] = true;
                        stack.push(nb);
                    }
                }
            }
            frag.sort_unstable();
            out.push(frag);
        }
        out
    }

    /// The symmetrized ring set, each ring as an ordered walk, in selection order.
    pub fn symmetrized(&self) -> Vec<Vec<usize>> {
        let mut degrees: Vec<i32> = self.adj.iter().map(|a| a.len() as i32).collect();
        let mut active = vec![true; self.bonds.len()];
        let mut invars: HashSet<Vec<usize>> = HashSet::new();
        let mut kept: Vec<Vec<usize>> = Vec::new();
        let mut extras: Vec<Vec<usize>> = Vec::new();

        for frag in self.fragments() {
            // the cyclomatic number of the fragment: how many rings it should yield
            let n_expected = (self.fragment_bond_count(&frag) + 1).saturating_sub(frag.len());
            if n_expected < 1 {
                continue;
            }
            let mut cands: Vec<Vec<usize>> = Vec::new();
            let mut changed: VecDeque<usize> =
                frag.iter().copied().filter(|&a| degrees[a] < 2).collect();
            let mut done = vec![false; self.n];
            let mut n_done = 0usize;

            while n_done + 3 <= frag.len() {
                while let Some(c) = changed.pop_front() {
                    if !done[c] {
                        done[c] = true;
                        n_done += 1;
                        self.trim_bonds(c, &mut changed, &mut degrees, &mut active);
                    }
                }
                let d2 = self.pick_d2_nodes(&frag, &degrees, &active);
                if !d2.is_empty() {
                    self.find_rings_d2_nodes(
                        &mut cands, &mut invars, &d2, &mut degrees, &mut active,
                    );
                    for &a in &d2 {
                        if !done[a] {
                            done[a] = true;
                            n_done += 1;
                        }
                        self.trim_bonds(a, &mut changed, &mut degrees, &mut active);
                    }
                } else if n_done + 3 <= frag.len() {
                    let d3 = frag.iter().copied().find(|&a| degrees[a] == 3);
                    match d3 {
                        None => break,
                        Some(c) => {
                            self.find_rings_d3_node(&mut cands, &mut invars, c, &active);
                            if !done[c] {
                                done[c] = true;
                                n_done += 1;
                            }
                            self.trim_bonds(c, &mut changed, &mut degrees, &mut active);
                        }
                    }
                } else {
                    break;
                }
            }

            // normalized before the cleanup, as `findSSSR` does, so extras carry it too
            for ring in &mut cands {
                normalize_ring(ring);
            }

            // the cleanup only runs when the fragment yielded more rings than its cyclomatic
            // number allows; otherwise the rings stay in the order they were found
            if cands.len() > n_expected {
                let (frag_kept, frag_extras) = self.remove_extra_rings(cands);
                kept.extend(frag_kept);
                extras.extend(frag_extras);
            } else {
                kept.extend(cands);
            }
        }

        let mut out = kept.clone();

        // an extra rejoins when it is the same size as a kept ring, shares a bond with it, and
        // supplies every bond that ring alone provides
        let kept_bonds: Vec<HashSet<usize>> =
            kept.iter().map(|r| self.ring_bond_set(r)).collect();
        let mut bond_counts: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for bs in &kept_bonds {
            for &b in bs {
                *bond_counts.entry(b).or_insert(0) += 1;
            }
        }
        for extra in &extras {
            let ebs = self.ring_bond_set(extra);
            for (ri, r) in kept.iter().enumerate() {
                if r.len() != extra.len() {
                    continue;
                }
                let mut share = false;
                let mut covers = true;
                for &b in &kept_bonds[ri] {
                    if ebs.contains(&b) {
                        share = true;
                    } else if bond_counts.get(&b) == Some(&1) {
                        covers = false;
                    }
                }
                if share && covers {
                    out.push(extra.clone());
                    break;
                }
            }
        }

        // Keep the order rings were selected in: kept rings in candidate order, then the
        // symmetrised extras. `max(key=len)` downstream takes the first maximum, so order decides
        // which ring is the pinned core when sizes tie.
        //
        // Each ring is emitted as the walk it was found as, not as a sorted atom list: the flat-ring
        // torsions read consecutive quadruples around a ring, so the traversal is part of the value.
        // Deduplication still keys on the sorted invariant, so a cycle found by two routes collapses.
        let mut seen: HashSet<Vec<usize>> = HashSet::new();
        let mut rings: Vec<Vec<usize>> = Vec::new();
        for r in &out {
            if seen.insert(invariant(r)) {
                rings.push(r.clone());
            }
        }
        rings
    }
}

/// Are `a` and `b` adjacent in the cyclic order of `ring`?
fn ring_adjacent(ring: &[usize], a: usize, b: usize) -> bool {
    let n = ring.len();
    for i in 0..n {
        let (x, y) = (ring[i], ring[(i + 1) % n]);
        if (x == a && y == b) || (x == b && y == a) {
            return true;
        }
    }
    false
}

/// Bonds lying on some cycle, by index into `bonds` — RDKit's `Bond::IsInRing`.
///
/// A bond is on a cycle exactly when it is not a cut edge. Degree-1 trimming is not enough: a
/// bridge joining two ring systems keeps both its atoms at degree 2 and survives trimming, but
/// lies on no cycle.
pub fn ring_bond_flags(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<bool> {
    if n_atoms == 0 || bonds.is_empty() {
        return vec![false; bonds.len()];
    }
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_atoms];
    for (bi, &(a, b)) in bonds.iter().enumerate() {
        adj[a].push((b, bi));
        adj[b].push((a, bi));
    }
    let mut disc = vec![usize::MAX; n_atoms];
    let mut low = vec![usize::MAX; n_atoms];
    let mut in_ring = vec![true; bonds.len()];
    let mut timer = 0usize;

    // iterative DFS: (atom, bond we entered by, next neighbour index to visit)
    for start in 0..n_atoms {
        if disc[start] != usize::MAX {
            continue;
        }
        let mut stack: Vec<(usize, usize, usize)> = vec![(start, usize::MAX, 0)];
        disc[start] = timer;
        low[start] = timer;
        timer += 1;
        while let Some(&mut (u, via, ref mut next)) = stack.last_mut() {
            if *next < adj[u].len() {
                let (v, bi) = adj[u][*next];
                *next += 1;
                if bi == via {
                    continue;
                }
                if disc[v] == usize::MAX {
                    disc[v] = timer;
                    low[v] = timer;
                    timer += 1;
                    stack.push((v, bi, 0));
                } else {
                    low[u] = low[u].min(disc[v]);
                }
            } else {
                stack.pop();
                if let Some(&mut (p, pv, _)) = stack.last_mut() {
                    low[p] = low[p].min(low[u]);
                    if low[u] > disc[p] {
                        in_ring[via] = false; // cut edge
                    }
                    let _ = pv;
                }
            }
        }
    }
    in_ring
}

/// Indices of the non-dative bonds. RDKit perceives rings with `includeDativeBonds=false`, so a
/// ligand→metal dative bond is invisible to SSSR: it is neither traversed to close a cycle nor
/// counted as a ring bond.
fn non_dative_indices(bonds: &[(usize, usize)], dative_donor: &[Option<usize>]) -> Vec<usize> {
    (0..bonds.len())
        .filter(|&bi| dative_donor.get(bi).copied().flatten().is_none())
        .collect()
}

/// [`symmetrized_sssr`] over the graph with dative bonds removed (RDKit's default
/// `includeDativeBonds=false`). Rings are atom lists, so dropping bonds keeps their indices valid.
pub fn symmetrized_sssr_masked(
    n_atoms: usize,
    bonds: &[(usize, usize)],
    dative_donor: &[Option<usize>],
) -> Vec<Vec<usize>> {
    let sub: Vec<(usize, usize)> = non_dative_indices(bonds, dative_donor)
        .into_iter()
        .map(|bi| bonds[bi])
        .collect();
    symmetrized_sssr(n_atoms, &sub)
}

/// [`ring_bond_flags`] with dative bonds removed from the graph, returning flags in the original
/// bond order — a dative bond is never a ring bond, and it cannot close a cycle for the others.
pub fn ring_bond_flags_masked(
    n_atoms: usize,
    bonds: &[(usize, usize)],
    dative_donor: &[Option<usize>],
) -> Vec<bool> {
    let kept = non_dative_indices(bonds, dative_donor);
    let sub: Vec<(usize, usize)> = kept.iter().map(|&bi| bonds[bi]).collect();
    let sub_flags = ring_bond_flags(n_atoms, &sub);
    let mut out = vec![false; bonds.len()];
    for (k, &bi) in kept.iter().enumerate() {
        out[bi] = sub_flags[k];
    }
    out
}

/// Convenience wrapper matching `rings::symmetrized_sssr`.
pub fn symmetrized_sssr(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
    if n_atoms == 0 || bonds.is_empty() {
        return Vec::new();
    }
    Rings::new(n_atoms, bonds).symmetrized()
}
