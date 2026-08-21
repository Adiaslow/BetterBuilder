//! Subgraph isomorphism in RDKit's traversal order.
//!
//! Matching a torsion pattern usually has several solutions, and `getExperimentalTorsions` keeps
//! the first one reported for a given central bond and locks the bond against later patterns. The
//! chosen quadruple is therefore decided by the order matches are discovered, not by any property
//! of the match, so the traversal order is part of the contract rather than an implementation
//! detail.
//!
//! RDKit vendors its own VF2 (`Code/GraphMol/Substruct/vf2.hpp`) and calls it with
//! `sortNodes = false`, so query nodes are taken in plain index order. Candidate target nodes come
//! from the neighbour list of an already-mapped atom, in adjacency order — which for RDKit's
//! `adjacency_list` is the order bonds were added. Node-frequency ordering, and the terminal-set
//! pruning behind `RDK_VF2_PRUNING`, are both compiled out.
//!
//! Matching is monomorphism, not induced subgraph isomorphism: every query bond must be present in
//! the target, but the target may carry bonds the query does not mention.

/// A graph as adjacency lists in edge-insertion order.
pub struct Graph {
    /// per node, its `(neighbour, edge index)` pairs in the order the edges were added
    adj: Vec<Vec<(usize, usize)>>,
}

impl Graph {
    /// Build from an edge list; `n` is the node count and edges are taken in index order.
    pub fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut adj = vec![Vec::new(); n];
        for (ei, &(a, b)) in edges.iter().enumerate() {
            adj[a].push((b, ei));
            adj[b].push((a, ei));
        }
        Graph { adj }
    }

    fn len(&self) -> usize {
        self.adj.len()
    }

    fn degree(&self, v: usize) -> usize {
        self.adj[v].len()
    }

    fn edge_between(&self, a: usize, b: usize) -> Option<usize> {
        self.adj[a].iter().find(|&&(n, _)| n == b).map(|&(_, e)| e)
    }
}

struct State<'a> {
    q: &'a Graph,
    t: &'a Graph,
    /// query node -> target node
    core_q: Vec<Option<usize>>,
    core_t: Vec<Option<usize>>,
    /// 0 when outside the terminal set, else the core length at which the node entered it
    term_q: Vec<usize>,
    term_t: Vec<usize>,
    core_len: usize,
    t_len_q: usize,
    t_len_t: usize,
}

impl<'a> State<'a> {
    fn new(q: &'a Graph, t: &'a Graph) -> Self {
        State {
            core_q: vec![None; q.len()],
            core_t: vec![None; t.len()],
            term_q: vec![0; q.len()],
            term_t: vec![0; t.len()],
            q,
            t,
            core_len: 0,
            t_len_q: 0,
            t_len_t: 0,
        }
    }

    /// The query node to extend with, and its candidate target nodes, in RDKit's order.
    ///
    /// Once anything is mapped, the query node is the lowest-indexed unmapped node adjacent to the
    /// mapping, and candidates are drawn only from the neighbours of one already-mapped target
    /// atom — the image of the first mapped neighbour in the query node's own adjacency order.
    fn candidates(&self) -> Option<(usize, Vec<usize>)> {
        if self.t_len_q > self.core_len && self.t_len_t > self.core_len {
            let n1 = (0..self.q.len())
                .find(|&i| self.core_q[i].is_none() && self.term_q[i] != 0)?;
            let anchor = self.q.adj[n1].iter().find_map(|&(u, _)| self.core_q[u])?;
            let cands = self.t.adj[anchor]
                .iter()
                .filter(|&&(v, _)| self.core_t[v].is_none())
                .map(|&(v, _)| v)
                .collect();
            Some((n1, cands))
        } else {
            let n1 = (0..self.q.len()).find(|&i| self.core_q[i].is_none())?;
            let cands = (0..self.t.len()).filter(|&v| self.core_t[v].is_none()).collect();
            Some((n1, cands))
        }
    }

    fn feasible<V, E>(&self, n1: usize, n2: usize, vc: &V, ec: &E) -> bool
    where
        V: Fn(usize, usize) -> bool,
        E: Fn(usize, usize) -> bool,
    {
        if self.q.degree(n1) > self.t.degree(n2) || !vc(n1, n2) {
            return false;
        }
        for &(other1, e1) in &self.q.adj[n1] {
            if let Some(other2) = self.core_q[other1] {
                match self.t.edge_between(n2, other2) {
                    Some(e2) if ec(e1, e2) => {}
                    _ => return false,
                }
            }
        }
        true
    }

    fn add_pair(&mut self, n1: usize, n2: usize) {
        self.core_len += 1;
        let len = self.core_len;
        for (term, len_acc, node, adj) in [
            (&mut self.term_q, &mut self.t_len_q, n1, &self.q.adj),
            (&mut self.term_t, &mut self.t_len_t, n2, &self.t.adj),
        ] {
            if term[node] == 0 {
                term[node] = len;
                *len_acc += 1;
            }
            for &(other, _) in &adj[node] {
                if term[other] == 0 {
                    term[other] = len;
                    *len_acc += 1;
                }
            }
        }
        self.core_q[n1] = Some(n2);
        self.core_t[n2] = Some(n1);
    }

    fn backtrack(&mut self, n1: usize, n2: usize) {
        let len = self.core_len;
        for (term, len_acc, node, adj) in [
            (&mut self.term_q, &mut self.t_len_q, n1, &self.q.adj),
            (&mut self.term_t, &mut self.t_len_t, n2, &self.t.adj),
        ] {
            if term[node] == len {
                term[node] = 0;
                *len_acc -= 1;
            }
            for &(other, _) in &adj[node] {
                if term[other] == len {
                    term[other] = 0;
                    *len_acc -= 1;
                }
            }
        }
        self.core_q[n1] = None;
        self.core_t[n2] = None;
        self.core_len -= 1;
    }

    /// Returns true once `limit` matches have been collected, unwinding the search.
    fn search<V, E>(&mut self, out: &mut Vec<Vec<usize>>, limit: usize, vc: &V, ec: &E) -> bool
    where
        V: Fn(usize, usize) -> bool,
        E: Fn(usize, usize) -> bool,
    {
        if self.core_len == self.q.len() {
            out.push(self.core_q.iter().map(|m| m.expect("core is full")).collect());
            return limit != 0 && out.len() >= limit;
        }
        if self.q.len() > self.t.len() || self.t_len_q > self.t_len_t {
            return false;
        }
        let Some((n1, cands)) = self.candidates() else {
            return false;
        };
        for n2 in cands {
            if self.feasible(n1, n2, vc, ec) {
                self.add_pair(n1, n2);
                if self.search(out, limit, vc, ec) {
                    return true;
                }
                self.backtrack(n1, n2);
            }
        }
        false
    }
}

/// Every mapping of `q` onto `t`, as query node -> target node, in discovery order.
///
/// Mappings are not deduplicated: `getExperimentalTorsions` matches with `uniquify = false`, and
/// symmetry-equivalent images are what the discovery order is there to arbitrate between. `limit`
/// caps the number returned, matching `SubstructMatchParams::maxMatches`; 0 means no cap.
pub fn all_matches<V, E>(q: &Graph, t: &Graph, vc: V, ec: E, limit: usize) -> Vec<Vec<usize>>
where
    V: Fn(usize, usize) -> bool,
    E: Fn(usize, usize) -> bool,
{
    let mut out = Vec::new();
    if q.len() == 0 || q.len() > t.len() {
        return out;
    }
    State::new(q, t).search(&mut out, limit, &vc, &ec);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A path query on a star: the centre is the only degree-3 atom, so candidate order is the
    /// centre's adjacency order.
    #[test]
    fn candidates_follow_adjacency_order() {
        let q = Graph::new(2, &[(0, 1)]);
        let t = Graph::new(4, &[(0, 1), (0, 2), (0, 3)]);
        let m = all_matches(&q, &t, |_, _| true, |_, _| true, 0);
        // query node 0 walks the target in index order; node 1 follows node 0's neighbours
        assert_eq!(m[0], vec![0, 1]);
        assert_eq!(m[1], vec![0, 2]);
        assert_eq!(m[2], vec![0, 3]);
        // then the leaves map to the centre
        assert_eq!(m[3], vec![1, 0]);
    }

    /// Reversing the edge insertion order reverses the order matches are found in.
    #[test]
    fn edge_order_decides_which_match_is_first() {
        let q = Graph::new(2, &[(0, 1)]);
        let a = all_matches(&Graph::new(2, &[(0, 1)]), &Graph::new(3, &[(0, 1), (0, 2)]), |_, _| true, |_, _| true, 1);
        let b = all_matches(&q, &Graph::new(3, &[(0, 2), (0, 1)]), |_, _| true, |_, _| true, 1);
        assert_eq!(a[0], vec![0, 1]);
        assert_eq!(b[0], vec![0, 2]);
    }

    /// Matching is monomorphism: an extra target bond between mapped atoms does not disqualify.
    #[test]
    fn extra_target_bonds_are_allowed() {
        let q = Graph::new(3, &[(0, 1), (1, 2)]);
        let t = Graph::new(3, &[(0, 1), (1, 2), (0, 2)]);
        assert!(!all_matches(&q, &t, |_, _| true, |_, _| true, 0).is_empty());
    }

    #[test]
    fn predicates_filter() {
        let q = Graph::new(1, &[]);
        let t = Graph::new(3, &[(0, 1), (1, 2)]);
        let m = all_matches(&q, &t, |_, ti| ti == 2, |_, _| true, 0);
        assert_eq!(m, vec![vec![2]]);
    }
}
