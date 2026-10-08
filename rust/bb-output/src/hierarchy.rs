//! The conformer-hierarchy builder — a faithful port of `mol2db2_py3_strain/hierarchy.py`.
//!
//! Produces the `X`/`R`/`C`/`S` coordinate hierarchy from the typed molecule + conformations:
//! independently-moving rigid structures ([`get_rigid_structures`]), per-atom position clustering and
//! conf/set assembly (`count_positions`, via [`crate::buckets`]), and the rigid core
//! (`find_rigid_heavy`). Every intermediate is exposed on [`Built`] so it can be gated against
//! the container's own `hierarchy.Hierarchy` internals, not just the final db2.

use crate::buckets::Buckets;
use crate::mol::TypedMol;
use crate::unionfind::UnionFind;
use std::collections::HashMap;

/// The builder's results: the intermediates (for gating) plus everything the db2 records need.
pub struct Built {
    pub rigid_structures: Vec<usize>,
    pub rigid_structure_count: usize,
    pub pos_count: Vec<usize>,
    pub num_confs: usize,
    pub out_atoms: usize,
    pub heavy_rigid_atom_nums: Vec<usize>,
    pub out_atom_orig_atom: Vec<usize>,
    pub out_atom_conf_num: Vec<usize>,
    pub out_atom_xyz: Vec<[f64; 3]>,
    pub conf_num_atom_list: Vec<[usize; 2]>,
    /// set index (input conf) -> 0-based conf numbers
    pub set_to_confs: Vec<Vec<usize>>,
    pub broken_sets: Vec<usize>,
}

fn is_heavy(sybyl: &str) -> bool {
    !sybyl.contains('H')
}

/// `_getRigidStructures`: map each atom to an independently-moving rigid structure. Atoms in a graph
/// cycle (ring system), joined by a non-single bond, or a degree-1 leaf with its neighbour, move
/// together; intersecting cycles merge into one structure.
pub fn get_rigid_structures(natoms: usize, atom_bonds: &[Vec<(usize, String)>]) -> (Vec<usize>, usize) {
    let mut cycle_count = 0usize;
    let mut cycles: Vec<Vec<usize>> = vec![Vec::new(); natoms];
    let mut parent: Vec<i64> = vec![0; natoms];
    let mut visited: Vec<u8> = vec![0; natoms];

    // recursive DFS cycle finder (macrocycles are small; matches Python's recursion order)
    fn find_cycles(
        atom: usize,
        prev: i64,
        atom_bonds: &[Vec<(usize, String)>],
        parent: &mut [i64],
        visited: &mut [u8],
        cycles: &mut [Vec<usize>],
        cycle_count: &mut usize,
    ) {
        if visited[atom] == 2 {
            return;
        }
        if visited[atom] == 1 {
            // back edge to a gray node: tag the parent chain from prev up to atom
            let mut curr = prev;
            while curr != atom as i64 {
                cycles[curr as usize].push(*cycle_count);
                curr = parent[curr as usize];
            }
            cycles[curr as usize].push(*cycle_count);
            *cycle_count += 1;
            return;
        }
        visited[atom] = 1;
        parent[atom] = prev;
        for (other, _bt) in &atom_bonds[atom] {
            if *other as i64 == prev {
                continue;
            }
            find_cycles(*other, atom as i64, atom_bonds, parent, visited, cycles, cycle_count);
        }
        visited[atom] = 2;
    }
    if natoms > 0 {
        find_cycles(0, -1, atom_bonds, &mut parent, &mut visited, &mut cycles, &mut cycle_count);
    }

    // add rigid bonds as length-2 cycles: degree-1 leaves, and non-single bonds not already shared
    for atom in 0..natoms {
        if atom_bonds[atom].len() == 1 {
            let nb = atom_bonds[atom][0].0;
            cycles[atom].push(cycle_count);
            cycles[nb].push(cycle_count);
            cycle_count += 1;
            continue;
        }
        for (other, bt) in atom_bonds[atom].clone() {
            if other < atom {
                continue;
            }
            let shared = cycles[atom].iter().any(|c| cycles[other].contains(c));
            if !shared && bt != "1" {
                cycles[atom].push(cycle_count);
                cycles[other].push(cycle_count);
                cycle_count += 1;
            }
        }
    }

    // cycles that share an atom intersect
    let mut intersections: Vec<Vec<usize>> = vec![Vec::new(); cycle_count];
    for atom_cycles in &cycles {
        for i in 0..atom_cycles.len() {
            for j in (i + 1)..atom_cycles.len() {
                let (a, b) = (atom_cycles[i], atom_cycles[j]);
                if !intersections[a].contains(&b) {
                    intersections[a].push(b);
                }
                if !intersections[b].contains(&a) {
                    intersections[b].push(a);
                }
            }
        }
    }

    // merge intersecting cycles into structures via DFS
    let mut rigid_map = vec![0usize; cycle_count];
    let mut cvisited = vec![false; cycle_count];
    let mut structure_count = 0usize;
    fn merge(cycle: usize, prev: i64, inter: &[Vec<usize>], rigid_map: &mut [usize], cvisited: &mut [bool], sc: usize) {
        if cvisited[cycle] {
            return;
        }
        cvisited[cycle] = true;
        rigid_map[cycle] = sc;
        for &nb in &inter[cycle] {
            if nb as i64 == prev {
                continue;
            }
            merge(nb, cycle as i64, inter, rigid_map, cvisited, sc);
        }
    }
    for cycle in 0..cycle_count {
        if !cvisited[cycle] {
            merge(cycle, -1, &intersections, &mut rigid_map, &mut cvisited, structure_count);
            structure_count += 1;
        }
    }

    // atom -> structure (first cycle's structure, or a fresh singleton structure)
    let mut rigid_structures = vec![0usize; natoms];
    for atom in 0..natoms {
        if !cycles[atom].is_empty() {
            rigid_structures[atom] = rigid_map[cycles[atom][0]];
        } else {
            rigid_structures[atom] = structure_count;
            structure_count += 1;
        }
    }
    (rigid_structures, structure_count)
}

/// `_countPositions` + the conf/set assembly. Returns the output-coordinate arrays and set map.
#[allow(clippy::type_complexity)]
fn count_positions(
    xyz_data: &[Vec<[f64; 3]>],
    rigid_structures: &[usize],
    tolerance: f64,
) -> (Vec<usize>, usize, usize, Vec<usize>, Vec<usize>, Vec<[f64; 3]>, Vec<[usize; 2]>, Vec<Vec<usize>>) {
    let nmol2s = xyz_data.len();
    let natoms = xyz_data[0].len();
    let bucketer = Buckets::new(tolerance);

    // per-atom clustering into conf_clusters (insertion-ordered via `order`)
    let mut clusters: HashMap<Vec<usize>, crate::buckets::ClusterVal> = HashMap::new();
    let mut order: Vec<Vec<usize>> = Vec::new();
    let mut pos_count = Vec::with_capacity(natoms);
    let mut conf_num = 0usize;
    #[allow(clippy::needless_range_loop)] // `atom` indexes the inner (per-conf) dimension
    for atom in 0..natoms {
        let col: Vec<[f64; 3]> = (0..nmol2s).map(|m| xyz_data[m][atom]).collect();
        let (npos, cn) = bucketer.bucket(&col, atom, conf_num, &mut clusters, &mut order);
        conf_num = cn;
        pos_count.push(npos);
    }

    // total output coords
    let pos_total: usize = clusters.values().map(|(_, atoms, _)| atoms.len()).sum();
    let mut out_atom_orig_atom = vec![0usize; pos_total];
    let mut out_atom_conf_num = vec![0usize; pos_total];
    let mut out_atom_xyz = vec![[0.0; 3]; pos_total];
    let mut conf_num_atom_list: Vec<[usize; 2]> = Vec::new();
    let mut set_to_confs: Vec<Vec<usize>> = vec![Vec::new(); nmol2s];

    // process clusters, rigid (longest key) first; Python's sort is stable on -len(key)
    let mut keys: Vec<Vec<usize>> = order.clone();
    keys.sort_by_key(|k| std::cmp::Reverse(k.len())); // stable: ties keep insertion order
    let mut global = 0usize;
    let mut conf_act = 0usize;
    for key in &keys {
        let (_, atoms, xyzlist) = &clusters[key];
        // stable sort of (atom, rigidStructure) by rigidStructure
        let mut idx: Vec<usize> = (0..atoms.len()).collect();
        idx.sort_by_key(|&i| rigid_structures[atoms[i]]); // stable
        let mut rs_prev: i64 = -1;
        let mut started = false;
        for (out_i, &i) in idx.iter().enumerate() {
            let rs = rigid_structures[atoms[i]] as i64;
            if started && rs_prev == rs {
                // same conf
            } else {
                if started {
                    conf_num_atom_list[conf_act - 1][1] = global - 1;
                }
                conf_num_atom_list.push([global, 0]);
                for &setno in key {
                    set_to_confs[setno].push(conf_act);
                }
                conf_act += 1;
            }
            // NB: the source indexes `atoms[i]`/`xyzlist[i]` by the *sorted* position `out_i`
            out_atom_orig_atom[global] = atoms[out_i];
            out_atom_conf_num[global] = conf_act;
            out_atom_xyz[global] = xyzlist[out_i];
            global += 1;
            rs_prev = rs;
            started = true;
        }
        conf_num_atom_list[conf_act - 1][1] = global - 1;
    }

    (pos_count, conf_act, global, out_atom_orig_atom, out_atom_conf_num, out_atom_xyz, conf_num_atom_list, set_to_confs)
}

/// `_findRigidComponent` + `_findRigidHeavy`: the largest bonded set of single-position atoms, its
/// heavy members in ascending order (matching CPython's set iteration for the contiguous case).
fn find_rigid_heavy(pos_count: &[usize], atom_bonds: &[Vec<(usize, String)>], atom_type: &[String]) -> Vec<usize> {
    let mut uf = UnionFind::new();
    for atom in 0..pos_count.len() {
        if pos_count[atom] == 1 {
            for (other, _) in &atom_bonds[atom] {
                if pos_count[*other] == 1 {
                    uf.union(atom, *other);
                }
            }
        }
    }
    let lists = uf.to_lists();
    let mut best: Vec<usize> = Vec::new();
    for l in lists {
        if l.len() > best.len() {
            best = l;
        }
    }
    // Her `heavyRigidAtomNums` = the heavy members of `set(rigidComponent)` in CPython SET-ITERATION
    // order (hierarchy.py 417/428: `atomsAssigned = set(rigidComponent)`, then `for a in atomsAssigned`).
    // That order is hash-table-layout dependent, NOT sorted — matching sorted only for small/contiguous
    // cores. So iterate a faithful CPython-set simulation over `best` (the maxCluster, in its build
    // order = the insertion order), then keep the heavy atoms in that order.
    cpython_set_iter_order(&best)
        .into_iter()
        .filter(|&a| is_heavy(&atom_type[a]))
        .collect()
}

/// The iteration order of a CPython `set` built by inserting `items` in order (keys are atom indices;
/// `hash(int) == int`). A faithful port of CPython `setobject.c` (`set_add_entry` +
/// `set_table_resize` + `set_insert_clean`), validated against real CPython 3.10 on 20k random int
/// lists (0 mismatches).
/// Constants are CPython's: LINEAR_PROBES=9, PERTURB_SHIFT=5, PySet_MINSIZE=8.
fn cpython_set_iter_order(items: &[usize]) -> Vec<usize> {
    const LINEAR_PROBES: usize = 9;
    const PERTURB_SHIFT: u32 = 5;
    const MINSIZE: usize = 8;

    // set_insert_clean: place `key` in the first empty slot (no dummies during a clean rebuild).
    fn insert_clean(table: &mut [Option<usize>], mask: usize, key: usize) {
        let mut perturb = key;
        let mut i = key & mask;
        loop {
            if table[i].is_none() {
                table[i] = Some(key);
                return;
            }
            if i + LINEAR_PROBES <= mask {
                for j in 1..=LINEAR_PROBES {
                    if table[i + j].is_none() {
                        table[i + j] = Some(key);
                        return;
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb)) & mask;
        }
    }

    let mut mask = MINSIZE - 1;
    let mut table: Vec<Option<usize>> = vec![None; MINSIZE];
    let (mut fill, mut used) = (0usize, 0usize);

    for &key in items {
        let mut i = key & mask;
        let mut perturb = key;
        let mut inserted = false;
        'probe: loop {
            match table[i] {
                None => {
                    table[i] = Some(key);
                    fill += 1;
                    used += 1;
                    inserted = true;
                    break 'probe;
                }
                Some(k) if k == key => break 'probe, // already present
                _ => {}
            }
            if i + LINEAR_PROBES <= mask {
                for j in 1..=LINEAR_PROBES {
                    match table[i + j] {
                        None => {
                            table[i + j] = Some(key);
                            fill += 1;
                            used += 1;
                            inserted = true;
                            break 'probe;
                        }
                        Some(k) if k == key => break 'probe,
                        _ => {}
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb)) & mask;
        }
        if inserted && fill * 5 >= mask * 3 {
            let minused = if used > 50000 { used * 2 } else { used * 4 };
            let mut newsize = MINSIZE;
            while newsize <= minused {
                newsize <<= 1;
            }
            let mut newtable: Vec<Option<usize>> = vec![None; newsize];
            let newmask = newsize - 1;
            for slot in table.iter().flatten() {
                insert_clean(&mut newtable, newmask, *slot);
            }
            table = newtable;
            mask = newmask;
            fill = used; // no dummies
        }
    }
    table.into_iter().flatten().collect()
}

/// Build the full hierarchy, including clash-broken sets ([`crate::clash`]).
pub fn build(mol2: &TypedMol, tolerance: f64) -> Built {
    let natoms = mol2.atom_num.len();
    let (rigid_structures, rigid_structure_count) = get_rigid_structures(natoms, &mol2.atom_bonds);
    let (pos_count, num_confs, out_atoms, out_atom_orig_atom, out_atom_conf_num, out_atom_xyz, conf_num_atom_list, set_to_confs) =
        count_positions(&mol2.atom_xyz, &rigid_structures, tolerance);
    let heavy_rigid_atom_nums = find_rigid_heavy(&pos_count, &mol2.atom_bonds, &mol2.atom_type);
    let broken_sets = crate::clash::broken_sets(mol2);
    Built {
        rigid_structures,
        rigid_structure_count,
        pos_count,
        num_confs,
        out_atoms,
        heavy_rigid_atom_nums,
        out_atom_orig_atom,
        out_atom_conf_num,
        out_atom_xyz,
        conf_num_atom_list,
        set_to_confs,
        broken_sets,
    }
}
