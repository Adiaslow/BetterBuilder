//! Ring perception: a minimum cycle basis plus the symmetry-equivalent rings RDKit stores.
//!
//! The target is the set `SANITIZE_SYMMRINGS` leaves on a molecule, because that is what every
//! downstream query reads — ring-size SMARTS such as `r{9-}`, ring-atom unions, largest-ring
//! selection. That set is a minimum cycle basis together with same-size cycles that could
//! substitute for a basis ring without uncovering any bond, so it is deliberately larger than a
//! basis and its rings are linearly dependent.
//!
//! The basis is built by Horton candidate generation and greedy GF(2) elimination rather than by
//! transcribing RDKit's degree-2/degree-3 decomposition: same object, stated guarantees. Selection
//! among equal-weight candidates is made deterministic by ordering on the sorted atom list, so the
//! result does not depend on traversal order.

/// A fixed-width bitset over bond indices, for GF(2) independence testing.
#[derive(Clone, PartialEq, Eq)]
struct BondSet(Vec<u64>);

impl BondSet {
    fn new(nbits: usize) -> Self {
        Self(vec![0; nbits.div_ceil(64)])
    }
    fn set(&mut self, i: usize) {
        self.0[i / 64] |= 1u64 << (i % 64);
    }
    fn get(&self, i: usize) -> bool {
        self.0[i / 64] >> (i % 64) & 1 == 1
    }
    fn xor_with(&mut self, other: &Self) {
        for (a, b) in self.0.iter_mut().zip(&other.0) {
            *a ^= *b;
        }
    }
    fn is_empty(&self) -> bool {
        self.0.iter().all(|w| *w == 0)
    }
    fn first_set(&self) -> Option<usize> {
        self.0
            .iter()
            .enumerate()
            .find(|(_, w)| **w != 0)
            .map(|(i, w)| i * 64 + w.trailing_zeros() as usize)
    }
    fn count(&self) -> usize {
        self.0.iter().map(|w| w.count_ones() as usize).sum()
    }
}

struct Graph {
    n: usize,
    /// `adj[a]` is `(neighbour, bond index)`.
    adj: Vec<Vec<(usize, usize)>>,
    bonds: Vec<(usize, usize)>,
}

impl Graph {
    fn new(n_atoms: usize, bonds: &[(usize, usize)]) -> Self {
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

    /// Bonds that lie on some cycle: everything left after iteratively removing degree-1 atoms.
    fn ring_bonds(&self) -> Vec<bool> {
        let mut degree: Vec<usize> = self.adj.iter().map(Vec::len).collect();
        let mut alive_atom = vec![true; self.n];
        let mut alive_bond = vec![true; self.bonds.len()];
        let mut queue: Vec<usize> = (0..self.n).filter(|&i| degree[i] <= 1).collect();
        while let Some(a) = queue.pop() {
            if !alive_atom[a] {
                continue;
            }
            alive_atom[a] = false;
            for &(nb, bi) in &self.adj[a] {
                if alive_bond[bi] && alive_atom[nb] {
                    alive_bond[bi] = false;
                    degree[nb] -= 1;
                    if degree[nb] == 1 {
                        queue.push(nb);
                    }
                }
            }
        }
        alive_bond
    }
}

/// Shortest-path tree from `root` over live bonds: parent atom and parent bond per atom.
fn bfs_tree(g: &Graph, root: usize, live: &[bool]) -> (Vec<i64>, Vec<usize>, Vec<i64>) {
    let mut dist = vec![-1i64; g.n];
    let mut parent = vec![usize::MAX; g.n];
    let mut parent_bond = vec![-1i64; g.n];
    let mut q = std::collections::VecDeque::new();
    dist[root] = 0;
    q.push_back(root);
    while let Some(a) = q.pop_front() {
        for &(nb, bi) in &g.adj[a] {
            if !live[bi] || dist[nb] >= 0 {
                continue;
            }
            dist[nb] = dist[a] + 1;
            parent[nb] = a;
            parent_bond[nb] = bi as i64;
            q.push_back(nb);
        }
    }
    (dist, parent, parent_bond)
}

/// Atoms on the tree path from `a` up to `root`, `a` first.
fn path_to_root(parent: &[usize], a: usize) -> Vec<usize> {
    let mut out = vec![a];
    let mut cur = a;
    while parent[cur] != usize::MAX {
        cur = parent[cur];
        out.push(cur);
    }
    out
}

/// A candidate cycle: its atoms and its bond set.
struct Candidate {
    atoms: Vec<usize>,
    bonds: BondSet,
}

/// Horton candidates: for each root `v` and each live edge `(x, y)`, the cycle
/// `SP(v, x) + (x, y) + SP(y, v)`, kept when the two tree paths meet only at `v`.
fn horton_candidates(g: &Graph, live: &[bool]) -> Vec<Candidate> {
    let nb = g.bonds.len();
    let mut seen: std::collections::HashSet<Vec<u64>> = std::collections::HashSet::new();
    let mut out = Vec::new();
    let in_ring: Vec<usize> = (0..g.n)
        .filter(|&a| g.adj[a].iter().any(|&(_, bi)| live[bi]))
        .collect();

    for &root in &in_ring {
        let (dist, parent, parent_bond) = bfs_tree(g, root, live);
        for (bi, &(x, y)) in g.bonds.iter().enumerate() {
            if !live[bi] || dist[x] < 0 || dist[y] < 0 {
                continue;
            }
            let px = path_to_root(&parent, x);
            let py = path_to_root(&parent, y);
            // a genuine simple cycle: the two paths share only the root
            let sx: std::collections::HashSet<usize> = px.iter().copied().collect();
            let shared: Vec<usize> = py.iter().copied().filter(|a| sx.contains(a)).collect();
            if shared.len() != 1 || shared[0] != root {
                continue;
            }
            let mut bset = BondSet::new(nb);
            bset.set(bi);
            for &a in px.iter().take(px.len().saturating_sub(1)) {
                if parent_bond[a] >= 0 {
                    bset.set(parent_bond[a] as usize);
                }
            }
            for &a in py.iter().take(py.len().saturating_sub(1)) {
                if parent_bond[a] >= 0 {
                    bset.set(parent_bond[a] as usize);
                }
            }
            let mut atoms: Vec<usize> = px.iter().copied().chain(py.iter().copied()).collect();
            atoms.sort_unstable();
            atoms.dedup();
            if atoms.len() != bset.count() {
                continue; // not a simple cycle
            }
            if seen.insert(bset.0.clone()) {
                out.push(Candidate { atoms, bonds: bset });
            }
        }
    }
    out
}

/// RDKit's symmetrized SSSR: a minimum cycle basis plus admissible same-size substitutes.
///
/// `bonds` are `(i, j)` atom index pairs. Returns rings as sorted atom-index lists, themselves
/// sorted, so the result is a canonical set.
pub fn symmetrized_sssr(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
    if n_atoms == 0 || bonds.is_empty() {
        return Vec::new();
    }
    let g = Graph::new(n_atoms, bonds);
    let live = g.ring_bonds();
    if !live.iter().any(|b| *b) {
        return Vec::new();
    }

    let mut cands = horton_candidates(&g, &live);
    // deterministic order: smallest first, then by atom list
    cands.sort_by(|a, b| {
        a.atoms
            .len()
            .cmp(&b.atoms.len())
            .then_with(|| a.atoms.cmp(&b.atoms))
    });

    // cycle rank over the ring system
    let live_bonds = live.iter().filter(|b| **b).count();
    let ring_atoms: std::collections::HashSet<usize> = bonds
        .iter()
        .enumerate()
        .filter(|(i, _)| live[*i])
        .flat_map(|(_, &(a, b))| [a, b])
        .collect();
    let components = count_components(&g, &live, &ring_atoms);
    let rank = live_bonds + components - ring_atoms.len();

    // greedy GF(2) elimination: keep a candidate when it is independent of those already kept
    let nb = g.bonds.len();
    let mut reduced: Vec<(usize, BondSet)> = Vec::new(); // (pivot bond, reduced row)
    let mut basis: Vec<usize> = Vec::new();
    let mut extras: Vec<usize> = Vec::new();
    for (ci, c) in cands.iter().enumerate() {
        if basis.len() == rank {
            extras.push(ci);
            continue;
        }
        let mut row = c.bonds.clone();
        for (pivot, r) in &reduced {
            if row.get(*pivot) {
                row.xor_with(r);
            }
        }
        if row.is_empty() {
            extras.push(ci);
        } else {
            let pivot = row.first_set().expect("non-empty row has a pivot");
            reduced.push((pivot, row));
            basis.push(ci);
        }
    }
    let _ = nb;

    // symmetrization: an extra joins the set when it is the same size as some basis ring, shares a
    // bond with it, and supplies every bond that ring alone provides.
    let mut bond_counts = vec![0usize; g.bonds.len()];
    for &bi in &basis {
        for k in 0..g.bonds.len() {
            if cands[bi].bonds.get(k) {
                bond_counts[k] += 1;
            }
        }
    }
    let mut chosen: Vec<usize> = basis.clone();
    for &ei in &extras {
        let e = &cands[ei];
        for &bi in &basis {
            let r = &cands[bi];
            if r.atoms.len() != e.atoms.len() {
                continue;
            }
            let mut share = false;
            let mut covers_unique = true;
            for k in 0..g.bonds.len() {
                if !r.bonds.get(k) {
                    continue;
                }
                if e.bonds.get(k) {
                    share = true;
                } else if bond_counts[k] == 1 {
                    covers_unique = false;
                }
            }
            if share && covers_unique {
                chosen.push(ei);
                break;
            }
        }
    }

    let mut rings: Vec<Vec<usize>> = chosen
        .into_iter()
        .map(|i| {
            let mut a = cands[i].atoms.clone();
            a.sort_unstable();
            a
        })
        .collect();
    rings.sort();
    rings.dedup();
    rings
}

fn count_components(
    g: &Graph,
    live: &[bool],
    ring_atoms: &std::collections::HashSet<usize>,
) -> usize {
    let mut seen = vec![false; g.n];
    let mut n = 0;
    for &start in ring_atoms {
        if seen[start] {
            continue;
        }
        n += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(a) = stack.pop() {
            for &(nb, bi) in &g.adj[a] {
                if live[bi] && !seen[nb] {
                    seen[nb] = true;
                    stack.push(nb);
                }
            }
        }
    }
    n
}
