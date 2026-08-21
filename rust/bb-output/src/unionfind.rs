//! Union-find with path compression and union-by-rank — a faithful port of `unionfind2.py`.
//!
//! `to_lists` must reproduce the Python dict iteration order (insertion order of first-seen items),
//! because `_findRigidComponent` picks the *first* largest cluster from that order. Items are `usize`
//! (atom or conf indices); an insertion-ordered `Vec` mirrors the Python `_parents` dict's key order.

use std::collections::HashMap;

#[derive(Default)]
pub struct UnionFind {
    parents: HashMap<usize, usize>,
    ranks: HashMap<usize, usize>,
    order: Vec<usize>, // first-seen order, == Python dict key order
}

impl UnionFind {
    pub fn new() -> Self {
        Self::default()
    }

    /// Find the root of `name`, inserting it as a singleton if unseen, with path compression.
    pub fn find(&mut self, name: usize) -> usize {
        if let std::collections::hash_map::Entry::Vacant(e) = self.parents.entry(name) {
            e.insert(name);
            self.ranks.insert(name, 0);
            self.order.push(name);
            return name;
        }
        // collect the path to the root
        let mut path = vec![name];
        let mut parent = self.parents[&name];
        while parent != *path.last().unwrap() {
            path.push(parent);
            parent = self.parents[&parent];
        }
        for &item in &path[..path.len() - 1] {
            self.parents.insert(item, parent);
        }
        parent
    }

    /// Union the sets of `name` and `other` by rank (ties raise `name`'s root's rank).
    pub fn union(&mut self, name: usize, other: usize) {
        let one = self.find(name);
        let two = self.find(other);
        if one == two {
            return;
        }
        if self.ranks[&one] < self.ranks[&two] {
            self.parents.insert(one, two);
        } else {
            self.parents.insert(two, one);
            if self.ranks[&one] == self.ranks[&two] {
                *self.ranks.get_mut(&one).unwrap() += 1;
            }
        }
    }

    /// Groups of unioned items, each in insertion order, groups in first-seen-root order — matching
    /// `unionfind2.toLists`.
    pub fn to_lists(&mut self) -> Vec<Vec<usize>> {
        let order = self.order.clone();
        for &item in &order {
            self.find(item); // compress everything to direct roots
        }
        let mut lists: Vec<Vec<usize>> = Vec::new();
        let mut root_pos: HashMap<usize, usize> = HashMap::new();
        for &item in &order {
            let par = self.parents[&item];
            match root_pos.get(&par) {
                Some(&pos) => lists[pos].push(item),
                None => {
                    root_pos.insert(par, lists.len());
                    lists.push(vec![item]);
                }
            }
        }
        lists
    }
}
