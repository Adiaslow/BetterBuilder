//! [`TypedMol`] — BetterBuilder's neutral in-memory typed molecule, the single source both output
//! artifacts (mol2 and db2) serialize from: atoms with SYBYL types, bonds, and one or more
//! conformations of coordinates. This is *not* a mol2 parser — production fills it from our own
//! pipeline, and the gate fills it from the container's authoritative parse (via JSON). It carries
//! just what the hierarchy algorithm and the typing rules read, plus `bonded_to`, a faithful port of
//! `mol2.py`'s `bondedToActual`/`bondedToActualAll` used by `sybyl2dock` and the color rules.

use std::collections::BTreeMap;

#[derive(Clone, Debug, Default)]
pub struct TypedMol {
    pub name: String,
    pub prot_name: String,
    pub smiles: String,
    pub longname: String,
    pub atom_num: Vec<usize>,   // 1-based ids as written
    pub atom_name: Vec<String>,
    pub atom_type: Vec<String>, // SYBYL type
    /// adjacency: `atom_bonds[i]` = `(other_0based, bond_type)`, in bond-insertion order
    pub atom_bonds: Vec<Vec<(usize, String)>>,
    pub bond_num: Vec<usize>,
    pub bond_start: Vec<usize>, // 1-based
    pub bond_end: Vec<usize>,   // 1-based
    pub bond_type: Vec<String>,
    /// `atom_xyz[conf][atom_index]` = (x, y, z)
    pub atom_xyz: Vec<Vec<[f64; 3]>>,
    pub input_total_strain: Vec<f64>,
    pub input_max_strain: Vec<f64>,
    pub input_hydrogens: Vec<i32>,
}

impl TypedMol {
    /// Build the bond adjacency (`atom_bonds`) from the bond lists, exactly as `mol2.py` accumulates
    /// it while reading `@<TRIPOS>BOND` (each bond appended to both endpoints, in file order).
    pub fn build_adjacency(&mut self) {
        self.atom_bonds = vec![Vec::new(); self.atom_num.len()];
        for i in 0..self.bond_num.len() {
            let (a, b) = (self.bond_start[i] - 1, self.bond_end[i] - 1);
            let bt = self.bond_type[i].clone();
            self.atom_bonds[a].push((b, bt.clone()));
            self.atom_bonds[b].push((a, bt));
        }
    }

    /// True iff `atom_num` (1-based) has, at *exactly* `bonds_away` bonds, some atom whose SYBYL type
    /// contains `first_name`. Reproduces `bondedToActual` over `bondedToActualAll`.
    pub fn bonded_to(&self, atom_num: usize, first_name: &str, bonds_away: usize, last_bond: Option<&str>) -> bool {
        let levels = self.bonded_all(atom_num - 1, bonds_away, last_bond);
        levels
            .get(&bonds_away)
            .into_iter()
            .flatten()
            .any(|&a| self.atom_type[a].contains(first_name))
    }

    /// `bondedToActualAll`: BFS levels from `actual` (0-based) out to `bonds_away`; each atom is placed
    /// at the first level it is reached and excluded from later ones. When `last_bond` is set, the
    /// final hop must use a bond whose type contains it.
    fn bonded_all(&self, actual: usize, bonds_away: usize, last_bond: Option<&str>) -> BTreeMap<usize, Vec<usize>> {
        let mut levels: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        levels.insert(0, vec![actual]);
        // An atom is placed at the first level it is reached; `visited` marks every atom already placed
        // (in a prior level or the level being built) so membership is O(1), not an O(n) scan per edge.
        let mut visited = vec![false; self.atom_bonds.len()];
        visited[actual] = true;
        let mut checked = 0usize;
        while checked < bonds_away {
            let cur = levels.get(&checked).cloned().unwrap_or_default();
            let mut next = Vec::new();
            for start in cur {
                for (other, btype) in &self.atom_bonds[start] {
                    if visited[*other] {
                        continue;
                    }
                    if let Some(lb) = last_bond {
                        if checked + 1 == bonds_away && lb != "*" && !btype.contains(lb) {
                            continue;
                        }
                    }
                    visited[*other] = true;
                    next.push(*other);
                }
            }
            levels.insert(checked + 1, next);
            checked += 1;
        }
        levels
    }
}
