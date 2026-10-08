//! The topological distance matrix: shortest bond-count path between every pair of atoms.
//!
//! Reproduces `RDKit::MolOps::getDistanceMat(mol)` with its defaults (`useBO=false`,
//! `useAtomWts=false`), which `setTopolBounds` calls once and hands to `setLowerBoundVDW` (the
//! 1-5/1-6 van der Waals scaling keys off distances of exactly 4 and 5 bonds) and to `set14`.
//!
//! RDKit initialises off-diagonal entries to `LOCAL_INF = 1e8` and diagonals to 0, gives every bond
//! weight 1, and runs Floyd–Warshall. With unit weights the all-pairs shortest path is
//! algorithm-independent, so a breadth-first search from each atom yields identical values in less
//! time on the sparse molecular graph; unreachable pairs (disconnected components) keep `LOCAL_INF`.
//! The result is flat `n*n` row-major, matching RDKit's layout so callers index `d[i*n + j]`.

/// RDKit's `LOCAL_INF`, the value left in place for a pair with no connecting path.
pub const LOCAL_INF: f64 = 1e8;

/// The `n*n` row-major topological distance matrix over the post-`AddHs` atoms.
pub fn distance_matrix(n_atoms: usize, bonds: &[(usize, usize)]) -> Vec<f64> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];
    for &(a, b) in bonds {
        adj[a].push(b);
        adj[b].push(a);
    }

    let mut d = vec![LOCAL_INF; n_atoms * n_atoms];
    let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    for src in 0..n_atoms {
        // BFS from src; distances are exact shortest bond counts, matching Floyd-Warshall's values
        d[src * n_atoms + src] = 0.0;
        queue.clear();
        queue.push_back(src);
        while let Some(u) = queue.pop_front() {
            let du = d[src * n_atoms + u];
            for &v in &adj[u] {
                if d[src * n_atoms + v] == LOCAL_INF {
                    d[src * n_atoms + v] = du + 1.0;
                    queue.push_back(v);
                }
            }
        }
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_chain_distances() {
        // 0-1-2-3: distances are |i-j| bonds
        let d = distance_matrix(4, &[(0, 1), (1, 2), (2, 3)]);
        let at = |i: usize, j: usize| d[i * 4 + j];
        assert_eq!(at(0, 0), 0.0);
        assert_eq!(at(0, 1), 1.0);
        assert_eq!(at(0, 3), 3.0);
        assert_eq!(at(1, 3), 2.0);
        // symmetric
        assert_eq!(at(3, 0), 3.0);
    }

    #[test]
    fn ring_takes_the_shorter_arc() {
        // a 4-ring 0-1-2-3-0: opposite atoms are 2 bonds apart either way
        let d = distance_matrix(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
        let at = |i: usize, j: usize| d[i * 4 + j];
        assert_eq!(at(0, 2), 2.0);
        assert_eq!(at(1, 3), 2.0);
        assert_eq!(at(0, 1), 1.0);
    }

    #[test]
    fn disconnected_pairs_stay_local_inf() {
        // two separate fragments 0-1 and 2-3 (as in a salt written with '.')
        let d = distance_matrix(4, &[(0, 1), (2, 3)]);
        let at = |i: usize, j: usize| d[i * 4 + j];
        assert_eq!(at(0, 1), 1.0);
        assert_eq!(at(0, 2), LOCAL_INF);
        assert_eq!(at(1, 3), LOCAL_INF);
    }
}
