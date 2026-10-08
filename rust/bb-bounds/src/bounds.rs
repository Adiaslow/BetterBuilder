//! The distance-geometry bounds matrix — the pure-Rust port of `setTopolBounds`.
//!
//! Built in RDKit's order: 1-2 bounds from UFF bond lengths ([`set12`]), 1-3 bounds from bond angles
//! ([`set13`]), then (later milestones) 1-4, 1-5, and finally van der Waals lower bounds for every
//! remaining pair ([`set_lower_bound_vdw`]). A [`ComputedData`] accumulator threads bond lengths,
//! bond angles, and visited flags between the stages exactly as RDKit's does.
//!
//! The matrix uses RDKit's `BoundsMatrix` layout so the flat result matches `getData()`
//! byte-for-byte: for `i < j` the upper bound lives at `[i*n+j]` and the lower bound at `[j*n+i]`;
//! the diagonal is 0.

use crate::{dist, uff};
use bb_perceive::hybrid::Hybridization;
use bb_perceive::smarts_match::{Perceived, Stereo};
use bb_perceive::valence;
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

/// Tolerance added around a 1-2 bond length (`DIST12_DELTA` in RDKit).
const DIST12_DELTA: f64 = 0.01;
/// Tolerance around a 1-3 distance (`DIST13_TOL`).
const DIST13_TOL: f64 = 0.04;
/// General distance tolerance (`GEN_DIST_TOL`), used for 1-4 bounds.
const GEN_DIST_TOL: f64 = 0.06;
/// Smallest ring size treated as a macrocycle for the 1-4 heuristics.
const MIN_MACROCYCLE_RING_SIZE: usize = 9;
/// The default upper bound; also the "unconstrained" sentinel (`MAX_UPPER`).
const MAX_UPPER: f64 = 1000.0;
/// Extra flex for larger heteroatoms in conjugated 5-rings (empirical, RDKit's `0.2`).
const EXTRA_SQUISH: f64 = 0.2;
/// van der Waals scaling for atoms four bonds apart (`VDW_SCALE_15`).
const VDW_SCALE_15: f64 = 0.7;
/// Lower bound for a hydrogen-bond donor H / acceptor pair (`H_BOND_LENGTH`).
const H_BOND_LENGTH: f64 = 1.8;

/// A distance-bounds matrix, flat `n*n` row-major in RDKit's `BoundsMatrix` layout.
pub struct BoundsMat {
    pub n: usize,
    pub data: Vec<f64>,
    /// Set when `check_and_set` sees a matrix `_checkAndSetBounds` would reject — an upper not
    /// greater than its lower, or a non-positive lower bound with no positive incumbent. RDKit
    /// throws a `CHECK_INVARIANT` there (reported as "bad lower bound"), aborting the build; we
    /// record it and let the caller reject the molecule after the matrix is finished, which is
    /// behaviourally identical since a rejected molecule's bounds are discarded either way.
    pub bad: bool,
}

impl BoundsMat {
    /// `initBoundsMat` with RDKit's defaults: upper triangle `MAX_UPPER`, lower triangle and
    /// diagonal 0.
    pub fn new(n: usize) -> Self {
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                data[i * n + j] = MAX_UPPER; // upper bound of pair (i<j); lower [j*n+i] stays 0
            }
        }
        BoundsMat { n, data, bad: false }
    }

    #[inline]
    fn upper(&self, i: usize, j: usize) -> f64 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.data[lo * self.n + hi]
    }
    #[inline]
    fn lower(&self, i: usize, j: usize) -> f64 {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.data[hi * self.n + lo]
    }
    #[inline]
    fn set_upper(&mut self, i: usize, j: usize, v: f64) {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.data[lo * self.n + hi] = v;
    }
    #[inline]
    fn set_lower(&mut self, i: usize, j: usize, v: f64) {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        self.data[hi * self.n + lo] = v;
    }

    /// `_checkAndSetBounds` (the `setIfBetter=false` path, the only one `setTopolBounds` uses):
    /// tighten the lower bound and loosen the upper only conservatively, never past existing values.
    // Each `if/else if` mirrors RDKit's two-branch structure verbatim (BoundsMatrixBuilder.cpp, its
    // own `// FIX this`): the branches differ in *condition* (bound still at the sentinel vs. the new
    // value tightening within range) but share the same action. clippy sees only the identical bodies.
    #[allow(clippy::if_same_then_else)]
    fn check_and_set(&mut self, i: usize, j: usize, lb: f64, ub: f64) {
        let clb = self.lower(i, j);
        let cub = self.upper(i, j);
        // RDKit's two `CHECK_INVARIANT`s in `_checkAndSetBounds`: the upper must exceed the lower,
        // and the lower (new or incumbent) must be positive. A violation makes RDKit throw and
        // reject the molecule ("bad lower bound"); we flag it and stop touching this pair.
        // `!(ub > lb)` rather than `ub <= lb` is deliberate: like `CHECK_INVARIANT`, it must also
        // reject a NaN bound (for which every comparison is false), so the negation is load-bearing.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(ub > lb) || !(lb > DIST12_DELTA || clb > DIST12_DELTA) {
            self.bad = true;
            return;
        }
        if clb <= DIST12_DELTA {
            self.set_lower(i, j, lb);
        } else if lb < clb && lb > DIST12_DELTA {
            self.set_lower(i, j, lb); // conservative bound setting
        }
        if cub >= MAX_UPPER {
            self.set_upper(i, j, ub);
        } else if ub > cub && ub < MAX_UPPER {
            self.set_upper(i, j, ub);
        }
    }
}

/// The cis/trans class of a recorded 1-4 path (`Path14Configuration::Path14Type`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Path14Type {
    Cis,
    Trans,
    Other,
}

/// Accumulated per-molecule geometry threaded between the bound-setting stages, mirroring RDKit's
/// `ComputedData`. Symmetric per-bond matrices (`bond_angles`, `bond_adj`) initialise to -1.
pub struct ComputedData {
    n: usize,
    nb: usize,
    pub bond_lengths: Vec<f64>,
    bond_angles: Vec<f64>,
    bond_adj: Vec<i64>,
    visited12: Vec<bool>,
    visited13: Vec<bool>,
    visited14: Vec<bool>,
    /// recorded 1-4 paths as `(bid1, bid2, bid3, type)` — consumed by set15
    pub paths14: Vec<(usize, usize, usize, Path14Type)>,
    /// cis / trans path keys `bid1*nb*nb + bid2*nb + bid3` (both orderings inserted) — for set15
    cis_paths: HashSet<u64>,
    trans_paths: HashSet<u64>,
    /// per pair, whether set15 has already written a bound — `set15Atoms`
    set15_atoms: Vec<bool>,
}

impl ComputedData {
    fn new(n: usize, nb: usize) -> Self {
        ComputedData {
            n,
            nb,
            bond_lengths: vec![0.0; nb],
            bond_angles: vec![-1.0; nb * nb],
            bond_adj: vec![-1; nb * nb],
            visited12: vec![false; n * n],
            visited13: vec![false; n * n],
            visited14: vec![false; n * n],
            paths14: Vec::new(),
            cis_paths: HashSet::new(),
            trans_paths: HashSet::new(),
            set15_atoms: vec![false; n * n],
        }
    }
    #[inline]
    fn path_key(&self, b1: usize, b2: usize, b3: usize) -> u64 {
        let nb = self.nb as u64;
        b1 as u64 * nb * nb + b2 as u64 * nb + b3 as u64
    }
    /// Whether a pair `(a, b)` is already a 1-2, 1-3 or 1-4 distance — `visitedBound(pid, DIST14)`.
    #[inline]
    fn visited_upto14(&self, a: usize, b: usize) -> bool {
        let p = self.pid(a, b);
        self.visited12[p] || self.visited13[p] || self.visited14[p]
    }
    #[inline]
    fn bond_adj_atom(&self, b1: usize, b2: usize) -> usize {
        self.bond_adj[b1 * self.nb + b2] as usize
    }
    #[inline]
    fn insert_cis(&mut self, b1: usize, b2: usize, b3: usize) {
        let nb = self.nb as u64;
        self.cis_paths
            .insert(b1 as u64 * nb * nb + b2 as u64 * nb + b3 as u64);
        self.cis_paths
            .insert(b3 as u64 * nb * nb + b2 as u64 * nb + b1 as u64);
    }
    #[inline]
    fn insert_trans(&mut self, b1: usize, b2: usize, b3: usize) {
        let nb = self.nb as u64;
        self.trans_paths
            .insert(b1 as u64 * nb * nb + b2 as u64 * nb + b3 as u64);
        self.trans_paths
            .insert(b3 as u64 * nb * nb + b2 as u64 * nb + b1 as u64);
    }
    /// Whether a pair `(a, b)` is already a 1-2 or 1-3 distance — `visitedBound(pid, DIST13)`.
    #[inline]
    fn visited_upto13(&self, a: usize, b: usize) -> bool {
        let p = self.pid(a, b);
        self.visited12[p] || self.visited13[p]
    }
    #[inline]
    fn pid(&self, a: usize, b: usize) -> usize {
        a.min(b) * self.n + a.max(b)
    }
    #[inline]
    fn bond_angle(&self, b1: usize, b2: usize) -> f64 {
        self.bond_angles[b1 * self.nb + b2]
    }
    #[inline]
    fn set_bond_angle(&mut self, b1: usize, b2: usize, v: f64) {
        self.bond_angles[b1 * self.nb + b2] = v;
        self.bond_angles[b2 * self.nb + b1] = v;
    }
    #[inline]
    fn set_bond_adj(&mut self, b1: usize, b2: usize, v: i64) {
        self.bond_adj[b1 * self.nb + b2] = v;
        self.bond_adj[b2 * self.nb + b1] = v;
    }
}

/// Per-atom bond structure in the order RDKit's `getAtomBonds` yields (bond-index order): each entry
/// is `(bond_index, other_atom)`.
fn atom_bonds(mol: &Perceived) -> Vec<Vec<(usize, usize)>> {
    let n = mol.atomic_numbers.len();
    let mut out = vec![Vec::new(); n];
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        out[a].push((bi, b));
        out[b].push((bi, a));
    }
    out
}

/// Bond index between two atoms, if directly bonded.
fn bond_lookup(mol: &Perceived) -> HashMap<(usize, usize), usize> {
    let mut m = HashMap::new();
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        m.insert((a.min(b), a.max(b)), bi);
    }
    m
}

/// Whether an atom lies on a ring of exactly `size` atoms — RDKit's `isAtomInRingOfSize`.
fn atom_in_ring_of_size(mol: &Perceived, atom: usize, size: usize) -> bool {
    mol.rings
        .iter()
        .any(|r| r.len() == size && r.contains(&atom))
}

/// Whether bond `(a, b)` lies on a ring of exactly `size` atoms — RDKit's `isBondInRingOfSize`.
fn bond_in_ring_of_size(mol: &Perceived, a: usize, b: usize, size: usize) -> bool {
    mol.rings.iter().filter(|r| r.len() == size).any(|r| {
        let m = r.len();
        (0..m).any(|i| {
            let (x, y) = (r[i], r[(i + 1) % m]);
            (x == a && y == b) || (x == b && y == a)
        })
    })
}

/// RDKit's `isLargerSP2Atom`: heavier than Al, SP2, and in at least one ring.
fn is_larger_sp2(mol: &Perceived, atom: usize) -> bool {
    mol.atomic_numbers[atom] > 13
        && mol.hybridization[atom] == Hybridization::Sp2
        && mol.rings.iter().any(|r| r.contains(&atom))
}

/// The 1-3 distance from two bond lengths and the angle between them (law of cosines).
fn compute13_dist(d1: f64, d2: f64, angle: f64) -> f64 {
    (d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * angle.cos()).sqrt()
}

/// RDKit's `_setRingAngle`: the interior angle at a ring atom, by hybridization and ring size.
fn set_ring_angle(hyb: Hybridization, ring_size: usize) -> f64 {
    if (hyb == Hybridization::Sp2 && ring_size <= 8) || ring_size == 3 || ring_size == 4 {
        PI * (1.0 - 2.0 / ring_size as f64)
    } else if hyb == Hybridization::Sp3 {
        if ring_size == 5 {
            104.0 * PI / 180.0
        } else {
            109.5 * PI / 180.0
        }
    } else if hyb == Hybridization::Sp3d {
        105.0 * PI / 180.0
    } else if hyb == Hybridization::Sp3d2 {
        90.0 * PI / 180.0
    } else {
        120.0 * PI / 180.0
    }
}

/// `_set13BoundsHelper`: set the 1-3 distance bound between `aid1` and `aid3` from the angle at the
/// central atom, widening the tolerance for larger SP2 atoms.
fn set13_helper(
    mol: &Perceived,
    aid1: usize,
    aid2: usize,
    aid3: usize,
    angle: f64,
    data: &ComputedData,
    bounds: &mut BoundsMat,
    blookup: &HashMap<(usize, usize), usize>,
) {
    let bid1 = blookup[&(aid1.min(aid2), aid1.max(aid2))];
    let bid2 = blookup[&(aid2.min(aid3), aid2.max(aid3))];
    let mut dl = compute13_dist(data.bond_lengths[bid1], data.bond_lengths[bid2], angle);
    let mut tol = DIST13_TOL;
    for a in [aid1, aid2, aid3] {
        if is_larger_sp2(mol, a) {
            tol *= 2.0;
        }
    }
    let du = dl + tol;
    dl -= tol;
    bounds.check_and_set(aid1, aid3, dl, du);
}

/// 1-2 bounds (see [`ComputedData`] for the accumulator). Reproduces `set12Bounds`.
pub fn set12(mol: &Perceived, bounds: &mut BoundsMat, data: &mut ComputedData) {
    let n = mol.atomic_numbers.len();
    let mut squish = vec![false; n];
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        if mol.bond_conjugated[bi]
            && (mol.atomic_numbers[a] > 10 || mol.atomic_numbers[b] > 10)
            && bond_in_ring_of_size(mol, a, b, 5)
        {
            squish[a] = true;
            squish[b] = true;
        }
    }

    let detailed = uff::bond_rest_lengths_detailed(mol);
    for (bi, &(a, b)) in mol.bonds.iter().enumerate() {
        let bl = detailed[bi].length;
        data.bond_lengths[bi] = bl;
        if detailed[bi].has_params {
            let extra = if squish[a] || squish[b] {
                EXTRA_SQUISH
            } else {
                0.0
            };
            bounds.set_upper(a, b, bl + extra + DIST12_DELTA);
            bounds.set_lower(a, b, bl - extra - DIST12_DELTA);
        } else {
            bounds.set_upper(a, b, 1.5 * bl);
            bounds.set_lower(a, b, 0.5 * bl);
        }
        let pid = data.pid(a, b);
        data.visited12[pid] = true;
    }
}

/// 1-3 bounds from bond angles. Reproduces `set13Bounds`: ring atoms first (angles from the
/// ring-size table), then the remaining 1-3 paths, distinguishing fused-ring atoms already visited
/// from plain non-ring atoms.
pub fn set13(mol: &Perceived, bounds: &mut BoundsMat, data: &mut ComputedData) {
    let n = mol.atomic_numbers.len();
    let nb = mol.bonds.len();
    let ab = atom_bonds(mol);
    let blookup = bond_lookup(mol);

    // rings sorted by size (stable, preserving our SSSR selection order among equal sizes)
    let mut rings: Vec<&Vec<usize>> = mol.rings.iter().collect();
    rings.sort_by_key(|r| r.len());

    let mut visited = vec![0i64; n];
    let mut angle_taken = vec![0.0f64; n];
    let mut done_paths = vec![false; nb * nb];

    // first pass: ring atoms, angles from the ring-size table
    for ring in &rings {
        let r_size = ring.len();
        let mut aid1 = ring[r_size - 1];
        for i in 0..r_size {
            let aid2 = ring[i];
            let aid3 = if i == r_size - 1 {
                ring[0]
            } else {
                ring[i + 1]
            };
            let bid1 = blookup[&(aid1.min(aid2), aid1.max(aid2))];
            let bid2 = blookup[&(aid2.min(aid3), aid2.max(aid3))];
            let id1 = nb * bid1 + bid2;
            let id2 = nb * bid2 + bid1;
            let pid = data.pid(aid1, aid3);
            if !done_paths[id1] && !done_paths[id2] {
                let angle = set_ring_angle(mol.hybridization[aid2], r_size);
                if !data.visited12[pid] {
                    set13_helper(mol, aid1, aid2, aid3, angle, data, bounds, &blookup);
                    data.visited13[pid] = true;
                }
                data.set_bond_angle(bid1, bid2, angle);
                data.set_bond_adj(bid1, bid2, aid2 as i64);
                visited[aid2] += 1;
                angle_taken[aid2] += angle;
                done_paths[id1] = true;
                done_paths[id2] = true;
            }
            aid1 = aid2;
        }
    }

    // second pass: every atom's remaining 1-3 paths
    for aid2 in 0..n {
        let deg = ab[aid2].len();
        let n13 = deg * (deg.saturating_sub(1)) / 2;
        if n13 == visited[aid2] as usize {
            continue;
        }
        let ahyb = mol.hybridization[aid2];
        // RDKit selects the ring/fused vs plain-non-ring branch ONCE per atom, from `visited`
        // BEFORE looping over bond pairs (`if (visited[aid2] >= 1)` guards the whole pair loop).
        // Since we increment `visited[aid2]` as we set each pair, we must capture the decision here
        // rather than re-reading `visited` per pair — otherwise a non-ring atom's first pair (which
        // makes visited 1) would flip every later pair into the ring branch.
        let ring_atom = visited[aid2] >= 1;
        // pairs of incident bonds in bond order: outer i, inner j < i (matching getAtomBonds order)
        for i in 0..ab[aid2].len() {
            let (bid1, aid1) = ab[aid2][i];
            for &(bid2, aid3) in &ab[aid2][..i] {
                if data.bond_angle(bid1, bid2) >= 0.0 {
                    continue; // this bond pair already has an angle
                }
                let pid = data.pid(aid1, aid3);
                if ring_atom {
                    // fused-ring / ring-atom-in-between cases
                    let angle = if ahyb == Hybridization::Sp2 {
                        (2.0 * PI - angle_taken[aid2]) / (n13 as f64 - visited[aid2] as f64)
                    } else if ahyb == Hybridization::Sp3 {
                        if atom_in_ring_of_size(mol, aid2, 3) {
                            116.0 * PI / 180.0
                        } else if atom_in_ring_of_size(mol, aid2, 4) {
                            112.0 * PI / 180.0
                        } else {
                            109.5 * PI / 180.0
                        }
                    } else if deg == 5 {
                        105.0 * PI / 180.0
                    } else if deg == 6 {
                        135.0 * PI / 180.0
                    } else {
                        120.0 * PI / 180.0
                    };
                    if !data.visited12[pid] {
                        set13_helper(mol, aid1, aid2, aid3, angle, data, bounds, &blookup);
                        data.visited13[pid] = true;
                    }
                    data.set_bond_angle(bid1, bid2, angle);
                    data.set_bond_adj(bid1, bid2, aid2 as i64);
                    angle_taken[aid2] += angle;
                    visited[aid2] += 1;
                } else {
                    // plain non-ring atom: angle purely from hybridization
                    let angle = match ahyb {
                        Hybridization::Sp => PI,
                        Hybridization::Sp2 => 2.0 * PI / 3.0,
                        Hybridization::Sp3 => 109.5 * PI / 180.0,
                        Hybridization::Sp3d => 105.0 * PI / 180.0,
                        Hybridization::Sp3d2 => 135.0 * PI / 180.0,
                        _ => 120.0 * PI / 180.0,
                    };
                    if !data.visited12[pid] {
                        if deg <= 4 {
                            set13_helper(mol, aid1, aid2, aid3, angle, data, bounds, &blookup);
                        } else {
                            // high-degree centre: 180 deg max, arbitrary min
                            let dmax = data.bond_lengths[bid1] + data.bond_lengths[bid2];
                            bounds.check_and_set(aid1, aid3, 1.0, dmax * 1.2);
                        }
                        data.visited13[pid] = true;
                    }
                    data.set_bond_angle(bid1, bid2, angle);
                    data.set_bond_adj(bid1, bid2, aid2 as i64);
                    angle_taken[aid2] += angle;
                    visited[aid2] += 1;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// set14Bounds: 1-4 (torsion) distance bounds.
// ---------------------------------------------------------------------------------------------

/// The 1-4 distance for a cis (0°) torsion, from the three bond lengths and two angles.
fn compute14_dist_cis(d1: f64, d2: f64, d3: f64, a12: f64, a23: f64) -> f64 {
    let dx = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy = d3 * a23.sin() - d1 * a12.sin();
    (dx * dx + dy * dy).sqrt()
}

/// The 1-4 distance for a trans (180°) torsion.
fn compute14_dist_trans(d1: f64, d2: f64, d3: f64, a12: f64, a23: f64) -> f64 {
    let dx = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy = d3 * a23.sin() + d1 * a12.sin();
    (dx * dx + dy * dy).sqrt()
}

/// The 1-4 distance for an arbitrary torsion angle (used for the S-S 90° case).
fn compute14_dist_3d(d1: f64, d2: f64, d3: f64, a12: f64, a23: f64, tor: f64) -> f64 {
    // atom 1 in the xy-plane; atom 4 built in the plane then rotated about x by the torsion
    let (p1x, p1y) = (d1 * a12.cos(), d1 * a12.sin());
    let p4x = d2 - d3 * a23.cos();
    let p4y0 = d3 * a23.sin();
    let (p4y, p4z) = (p4y0 * tor.cos(), p4y0 * tor.sin());
    let (dx, dy, dz) = (p4x - p1x, p4y - p1y, p4z);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn stereo_is_definite(s: Stereo) -> bool {
    matches!(s, Stereo::Z | Stereo::E | Stereo::Cis | Stereo::Trans)
}
fn stereo_prefers_cis(s: Stereo) -> bool {
    matches!(s, Stereo::Z | Stereo::Cis)
}

/// Bundle of precomputed molecule data for the 1-4 pass, mirroring the accessors `set14` uses.
struct S14<'a> {
    mol: &'a Perceived,
    n: usize,
    nb: usize,
    ab: Vec<Vec<(usize, usize)>>, // per atom: (bond_index, other_atom), in bond-index order
    blookup: HashMap<(usize, usize), usize>,
    num_hs: Vec<usize>,
    num_bond_rings: Vec<usize>,
    bond_rings: Vec<Vec<usize>>,
    dmat: &'a [f64],
    force_trans_amides: bool,
    macrocycle14: bool,
}

impl<'a> S14<'a> {
    fn new(
        mol: &'a Perceived,
        dmat: &'a [f64],
        macrocycle14: bool,
        force_trans_amides: bool,
    ) -> Self {
        let n = mol.atomic_numbers.len();
        let nb = mol.bonds.len();
        let ab = atom_bonds(mol);
        let blookup = bond_lookup(mol);
        let num_hs: Vec<usize> = (0..n)
            .map(|a| {
                ab[a]
                    .iter()
                    .filter(|&&(_, o)| mol.atomic_numbers[o] == 1)
                    .count()
            })
            .collect();
        let bond_rings: Vec<Vec<usize>> = mol
            .rings
            .iter()
            .map(|ring| {
                let m = ring.len();
                (0..m)
                    .map(|i| {
                        let (x, y) = (ring[i], ring[(i + 1) % m]);
                        blookup[&(x.min(y), x.max(y))]
                    })
                    .collect()
            })
            .collect();
        let mut num_bond_rings = vec![0usize; nb];
        for br in &bond_rings {
            for &b in br {
                num_bond_rings[b] += 1;
            }
        }
        S14 {
            mol,
            n,
            nb,
            ab,
            blookup,
            num_hs,
            num_bond_rings,
            bond_rings,
            dmat,
            force_trans_amides,
            macrocycle14,
        }
    }

    #[inline]
    fn other_atom(&self, bid: usize, atom: usize) -> usize {
        let (x, y) = self.mol.bonds[bid];
        if x == atom {
            y
        } else {
            x
        }
    }
    #[inline]
    fn degree(&self, atom: usize) -> usize {
        self.ab[atom].len()
    }
    #[inline]
    fn bond_between(&self, a: usize, b: usize) -> bool {
        self.blookup.contains_key(&(a.min(b), a.max(b)))
    }
    #[inline]
    fn z(&self, atom: usize) -> u8 {
        self.mol.atomic_numbers[atom]
    }
    #[inline]
    fn hyb(&self, atom: usize) -> Hybridization {
        self.mol.hybridization[atom]
    }
    #[inline]
    fn bond_order(&self, bid: usize) -> u8 {
        self.mol.bond_order[bid]
    }
    #[inline]
    fn dmat_dist(&self, a: usize, b: usize) -> f64 {
        self.dmat[a.max(b) * self.n + a.min(b)]
    }

    /// `_getAtomStereo`: the bond stereo as seen from the specific 1-4 outer atoms.
    fn atom_stereo(&self, bid: usize, aid1: usize, aid4: usize) -> Stereo {
        let mut st = self.mol.bond_stereo[bid];
        let sa = self.mol.stereo_atoms[bid];
        if stereo_is_definite(st)
            && sa[0] >= 0
            && sa[1] >= 0
            && ((sa[0] != aid1 as i64) ^ (sa[1] != aid4 as i64))
        {
            st = match st {
                Stereo::Z => Stereo::E,
                Stereo::E => Stereo::Z,
                Stereo::Cis => Stereo::Trans,
                Stereo::Trans => Stereo::Cis,
                other => other,
            };
        }
        st
    }

    // --- amide/ester predicates ---

    fn is_carbonyl(&self, at: usize) -> bool {
        if self.z(at) == 6 && self.degree(at) > 2 {
            for &(bid, nbr) in &self.ab[at] {
                let an = self.z(nbr);
                if (an == 8 || an == 7) && self.bond_order(bid) == 2 {
                    return true;
                }
            }
        }
        false
    }

    /// `_checkAmideEster14`: 1(?)-2(N/O)-3(C=)-4(O/N), bnd1 single, bnd3 double.
    fn check_amide_ester_14(
        &self,
        bid1: usize,
        bid3: usize,
        a2: usize,
        a3: usize,
        a4: usize,
    ) -> bool {
        self.z(a3) == 6
            && self.bond_order(bid3) == 2
            && (self.z(a4) == 8 || self.z(a4) == 7)
            && self.bond_order(bid1) == 1
            && (self.z(a2) == 8 || (self.z(a2) == 7 && self.num_hs[a2] == 1))
    }

    /// `_checkAmideEster15`: a2 is O or NH1, bnd1 single, a3 is a carbonyl C, bnd3 single.
    fn check_amide_ester_15(&self, bid1: usize, bid3: usize, a2: usize, a3: usize) -> bool {
        let a2num = self.z(a2);
        (a2num == 8 || (a2num == 7 && self.num_hs[a2] == 1))
            && self.bond_order(bid1) == 1
            && self.z(a3) == 6
            && self.bond_order(bid3) == 1
            && self.is_carbonyl(a3)
    }

    /// `_checkMacrocycleAllInSameRingAmideEster14`.
    fn check_macro_all_amide_ester_14(&self, a1: usize, a2: usize, a3: usize, a4: usize) -> bool {
        if self.z(a3) != 6 {
            return false;
        }
        if (self.z(a2) == 7 || self.z(a2) == 8) && self.degree(a2) == 3 && self.degree(a3) == 3 {
            // atm2's third neighbour (not a1/a3) must be C or H via a single bond
            for &(bid, nbr) in &self.ab[a2] {
                if nbr != a1 && nbr != a3 {
                    if (self.z(nbr) != 6 && self.z(nbr) != 1) || self.bond_order(bid) != 1 {
                        return false;
                    }
                    break;
                }
            }
            // atm3's third neighbour (not a2/a4) must be a carbonyl O via a double bond
            for &(bid, nbr) in &self.ab[a3] {
                if nbr != a2 && nbr != a4 {
                    if self.z(nbr) != 8 || self.bond_order(bid) != 2 {
                        return false;
                    }
                    break;
                }
            }
            return true;
        }
        false
    }

    /// `_checkMacrocycleTwoInSameRingAmideEster14`.
    fn check_macro_two_amide_ester_14(
        &self,
        bid1: usize,
        bid3: usize,
        a1: usize,
        a2: usize,
        a3: usize,
        a4: usize,
    ) -> bool {
        self.z(a1) != 1
            && self.z(a3) == 6
            && self.bond_order(bid3) == 2
            && (self.z(a4) == 8 || self.z(a4) == 7)
            && self.bond_order(bid1) == 1
            && (self.z(a2) == 8 || self.z(a2) == 7)
    }

    // --- the six 1-4 setters ---

    fn record_14_path(&self, bid1: usize, bid2: usize, bid3: usize, data: &mut ComputedData) {
        let atm2 = data.bond_adj_atom(bid1, bid2);
        let atm3 = data.bond_adj_atom(bid2, bid3);
        if self.hyb(atm2) == Hybridization::Sp2 && self.hyb(atm3) == Hybridization::Sp2 {
            data.paths14.push((bid1, bid2, bid3, Path14Type::Cis));
            data.insert_cis(bid1, bid2, bid3);
        } else {
            data.paths14.push((bid1, bid2, bid3, Path14Type::Other));
        }
    }

    fn set_in_ring_14(
        &self,
        bid1: usize,
        bid2: usize,
        bid3: usize,
        ring_size: i64,
        bounds: &mut BoundsMat,
        data: &mut ComputedData,
    ) {
        let atm2 = data.bond_adj_atom(bid1, bid2);
        let atm3 = data.bond_adj_atom(bid2, bid3);
        let ahyb2 = self.hyb(atm2);
        let ahyb3 = self.hyb(atm3);
        let aid1 = self.other_atom(bid1, atm2);
        let aid4 = self.other_atom(bid3, atm3);
        if data.visited_upto13(aid1, aid4) {
            return;
        }
        if self.dmat_dist(aid1, aid4) < 2.9 {
            return;
        }
        let (bl1, bl2, bl3) = (
            data.bond_lengths[bid1],
            data.bond_lengths[bid2],
            data.bond_lengths[bid3],
        );
        let ba12 = data.bond_angle(bid1, bid2);
        let ba23 = data.bond_angle(bid2, bid3);
        let stype = self.atom_stereo(bid2, aid1, aid4);
        let mut prefer_cis = false;
        let mut prefer_trans = false;
        if ring_size <= 8
            && ahyb2 == Hybridization::Sp2
            && ahyb3 == Hybridization::Sp2
            && stype != Stereo::E
            && stype != Stereo::Trans
        {
            if self.num_bond_rings[bid2] > 1 {
                if self.num_bond_rings[bid1] == 1 && self.num_bond_rings[bid3] == 1 {
                    for br in &self.bond_rings {
                        if br.contains(&bid1) {
                            if br.contains(&bid3) {
                                prefer_cis = true;
                            }
                            break;
                        }
                    }
                }
            } else {
                prefer_cis = true;
            }
        } else if stype == Stereo::Z || stype == Stereo::Cis {
            prefer_cis = true;
        } else if stype == Stereo::E || stype == Stereo::Trans {
            prefer_trans = true;
        }

        let ptype;
        let (mut dl, mut du);
        if prefer_cis {
            ptype = Path14Type::Cis;
            data.insert_cis(bid1, bid2, bid3);
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
            du = dl + 2.0 * GEN_DIST_TOL;
        } else if prefer_trans {
            ptype = Path14Type::Trans;
            data.insert_trans(bid1, bid2, bid3);
            dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
            du = dl + 2.0 * GEN_DIST_TOL;
        } else {
            ptype = Path14Type::Other;
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
            du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
            if du < dl {
                std::mem::swap(&mut du, &mut dl);
            }
            if (du - dl).abs() < DIST12_DELTA {
                dl -= GEN_DIST_TOL;
                du += GEN_DIST_TOL;
            }
        }
        data.paths14.push((bid1, bid2, bid3, ptype));
        let p = data.pid(aid1, aid4);
        data.visited14[p] = true;
        bounds.check_and_set(aid1, aid4, dl, du);
    }

    fn set_two_in_same_ring_14(
        &self,
        bid1: usize,
        bid2: usize,
        bid3: usize,
        bounds: &mut BoundsMat,
        data: &mut ComputedData,
    ) {
        let atm2 = data.bond_adj_atom(bid1, bid2);
        let atm3 = data.bond_adj_atom(bid2, bid3);
        let aid1 = self.other_atom(bid1, atm2);
        let aid4 = self.other_atom(bid3, atm3);
        if data.visited_upto13(aid1, aid4) {
            return;
        }
        if self.dmat_dist(aid1, aid4) < 2.9 {
            return;
        }
        // fused-ring guard (sf.net bug 2835784)
        if self.bond_between(aid1, atm3) || self.bond_between(aid4, atm2) {
            return;
        }
        let (bl1, bl2, bl3) = (
            data.bond_lengths[bid1],
            data.bond_lengths[bid2],
            data.bond_lengths[bid3],
        );
        let ba12 = data.bond_angle(bid1, bid2);
        let ba23 = data.bond_angle(bid2, bid3);
        let (mut dl, mut du);
        let ptype;
        if self.hyb(atm2) == Hybridization::Sp2 && self.hyb(atm3) == Hybridization::Sp2 {
            dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
            du = dl;
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
            ptype = Path14Type::Trans;
            data.insert_trans(bid1, bid2, bid3);
        } else {
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
            du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
            if du < dl {
                std::mem::swap(&mut dl, &mut du);
            }
            if (du - dl).abs() < DIST12_DELTA {
                dl -= GEN_DIST_TOL;
                du += GEN_DIST_TOL;
            }
            ptype = Path14Type::Other;
        }
        bounds.check_and_set(aid1, aid4, dl, du);
        data.paths14.push((bid1, bid2, bid3, ptype));
        let p = data.pid(aid1, aid4);
        data.visited14[p] = true;
    }

    fn set_macrocycle_two_in_same_ring_14(
        &self,
        bid1: usize,
        bid2: usize,
        bid3: usize,
        bounds: &mut BoundsMat,
        data: &mut ComputedData,
    ) {
        let atm2 = data.bond_adj_atom(bid1, bid2);
        let atm3 = data.bond_adj_atom(bid2, bid3);
        let aid1 = self.other_atom(bid1, atm2);
        let aid4 = self.other_atom(bid3, atm3);
        if data.visited_upto13(aid1, aid4) {
            return;
        }
        if self.dmat_dist(aid1, aid4) < 2.9 {
            return;
        }
        if self.bond_between(aid1, atm3) || self.bond_between(aid4, atm2) {
            return;
        }
        let (bl1, bl2, bl3) = (
            data.bond_lengths[bid1],
            data.bond_lengths[bid2],
            data.bond_lengths[bid3],
        );
        let ba12 = data.bond_angle(bid1, bid2);
        let ba23 = data.bond_angle(bid2, bid3);
        let (mut dl, mut du);
        let ptype;
        if self.check_macro_two_amide_ester_14(bid1, bid3, aid1, atm2, atm3, aid4)
            || self.check_macro_two_amide_ester_14(bid3, bid1, aid4, atm3, atm2, aid1)
        {
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
            ptype = Path14Type::Cis;
            data.insert_cis(bid1, bid2, bid3);
            du = dl;
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        } else {
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
            du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
            if du < dl {
                std::mem::swap(&mut dl, &mut du);
            }
            if (du - dl).abs() < DIST12_DELTA {
                dl -= GEN_DIST_TOL;
                du += GEN_DIST_TOL;
            }
            ptype = Path14Type::Other;
        }
        bounds.check_and_set(aid1, aid4, dl, du);
        data.paths14.push((bid1, bid2, bid3, ptype));
        let p = data.pid(aid1, aid4);
        data.visited14[p] = true;
    }

    // dl/du/ptype are assigned in every match arm; the initialisers mirror the C++ `double dl=0.0`
    // and are overwritten before use.
    #[allow(unused_assignments)]
    fn set_chain_14(
        &self,
        bid1: usize,
        bid2: usize,
        bid3: usize,
        bounds: &mut BoundsMat,
        data: &mut ComputedData,
    ) {
        let atm2 = data.bond_adj_atom(bid1, bid2);
        let atm3 = data.bond_adj_atom(bid2, bid3);
        let aid1 = self.other_atom(bid1, atm2);
        let aid4 = self.other_atom(bid3, atm3);
        if data.visited_upto13(aid1, aid4) {
            return;
        }
        let (bl1, bl2, bl3) = (
            data.bond_lengths[bid1],
            data.bond_lengths[bid2],
            data.bond_lengths[bid3],
        );
        let ba12 = data.bond_angle(bid1, bid2);
        let ba23 = data.bond_angle(bid2, bid3);
        let mut dl = 0.0;
        let mut du = 0.0;
        let mut ptype = Path14Type::Other;

        match self.bond_order(bid2) {
            2 => {
                if self.bond_order(bid1) == 2 || self.bond_order(bid3) == 2 {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                    du = dl + 2.0 * GEN_DIST_TOL;
                    ptype = Path14Type::Cis;
                    data.insert_cis(bid1, bid2, bid3);
                } else if stereo_is_definite(self.mol.bond_stereo[bid2]) {
                    let stype = self.atom_stereo(bid2, aid1, aid4);
                    if stereo_prefers_cis(stype) {
                        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                        du = dl + 2.0 * GEN_DIST_TOL;
                        ptype = Path14Type::Cis;
                        data.insert_cis(bid1, bid2, bid3);
                    } else {
                        du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        dl = du - GEN_DIST_TOL;
                        du += GEN_DIST_TOL;
                        ptype = Path14Type::Trans;
                        data.insert_trans(bid1, bid2, bid3);
                    }
                } else {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    if (du - dl).abs() < DIST12_DELTA {
                        dl -= GEN_DIST_TOL;
                        du += GEN_DIST_TOL;
                    }
                    ptype = Path14Type::Other;
                }
            }
            1 => {
                if self.z(atm2) == 16
                    && self.z(atm3) == 16
                    && self.degree(atm2) == 2
                    && self.degree(atm3) == 2
                {
                    dl = compute14_dist_3d(bl1, bl2, bl3, ba12, ba23, PI / 2.0) - GEN_DIST_TOL;
                    du = dl + 2.0 * GEN_DIST_TOL;
                    ptype = Path14Type::Other;
                } else if self.check_amide_ester_14(bid1, bid3, atm2, atm3, aid4)
                    || self.check_amide_ester_14(bid3, bid1, atm3, atm2, aid1)
                {
                    if self.force_trans_amides {
                        if (self.z(aid1) == 1
                            && self.z(atm2) == 7
                            && self.degree(atm2) == 3
                            && self.num_hs[atm2] == 1)
                            || (self.z(aid4) == 1
                                && self.z(atm3) == 7
                                && self.degree(atm3) == 3
                                && self.num_hs[atm3] == 1)
                        {
                            dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                            ptype = Path14Type::Trans;
                            data.insert_trans(bid1, bid2, bid3);
                        } else {
                            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                            ptype = Path14Type::Cis;
                            data.insert_cis(bid1, bid2, bid3);
                        }
                        du = dl;
                        dl -= GEN_DIST_TOL;
                        du += GEN_DIST_TOL;
                    } else {
                        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        ptype = Path14Type::Other;
                    }
                } else if self.check_amide_ester_15(bid1, bid3, atm2, atm3)
                    || self.check_amide_ester_15(bid3, bid1, atm3, atm2)
                {
                    if self.force_trans_amides {
                        if (self.z(aid1) == 1
                            && self.z(atm2) == 7
                            && self.degree(atm2) == 3
                            && self.num_hs[atm2] == 1)
                            || (self.z(aid4) == 1
                                && self.z(atm3) == 7
                                && self.degree(atm3) == 3
                                && self.num_hs[atm3] == 1)
                        {
                            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                            ptype = Path14Type::Cis;
                            data.insert_cis(bid1, bid2, bid3);
                        } else {
                            dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                            ptype = Path14Type::Trans;
                            data.insert_trans(bid1, bid2, bid3);
                        }
                        du = dl;
                        dl -= GEN_DIST_TOL;
                        du += GEN_DIST_TOL;
                    } else {
                        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        ptype = Path14Type::Other;
                    }
                } else {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    ptype = Path14Type::Other;
                }
            }
            _ => {
                dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                ptype = Path14Type::Other;
            }
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        bounds.check_and_set(aid1, aid4, dl, du);
        data.paths14.push((bid1, bid2, bid3, ptype));
        let p = data.pid(aid1, aid4);
        data.visited14[p] = true;
    }

    // dl/du/ptype are assigned in every reachable path; the initialisers mirror the C++ source.
    #[allow(unused_assignments)]
    fn set_macrocycle_all_in_same_ring_14(
        &self,
        bid1: usize,
        bid2: usize,
        bid3: usize,
        bounds: &mut BoundsMat,
        data: &mut ComputedData,
    ) {
        let atm2 = data.bond_adj_atom(bid1, bid2);
        let atm3 = data.bond_adj_atom(bid2, bid3);
        let aid1 = self.other_atom(bid1, atm2);
        let aid4 = self.other_atom(bid3, atm3);
        if data.visited_upto13(aid1, aid4) {
            return;
        }
        let (bl1, bl2, bl3) = (
            data.bond_lengths[bid1],
            data.bond_lengths[bid2],
            data.bond_lengths[bid3],
        );
        let ba12 = data.bond_angle(bid1, bid2);
        let ba23 = data.bond_angle(bid2, bid3);
        let mut dl = 0.0;
        let mut du = 0.0;
        let mut ptype = Path14Type::Other;
        let mut set_the_bound = true;

        match self.bond_order(bid2) {
            2 => {
                if self.bond_order(bid1) == 2 || self.bond_order(bid3) == 2 {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                    du = dl + 2.0 * GEN_DIST_TOL;
                    ptype = Path14Type::Cis;
                    data.insert_cis(bid1, bid2, bid3);
                } else if stereo_is_definite(self.mol.bond_stereo[bid2]) {
                    let stype = self.atom_stereo(bid2, aid1, aid4);
                    if stereo_prefers_cis(stype) {
                        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                        du = dl + 2.0 * GEN_DIST_TOL;
                        ptype = Path14Type::Cis;
                        data.insert_cis(bid1, bid2, bid3);
                    } else {
                        du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        dl = du - GEN_DIST_TOL;
                        du += GEN_DIST_TOL;
                        ptype = Path14Type::Trans;
                        data.insert_trans(bid1, bid2, bid3);
                    }
                } else {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    if (du - dl).abs() < DIST12_DELTA {
                        dl -= GEN_DIST_TOL;
                        du += GEN_DIST_TOL;
                    }
                    ptype = Path14Type::Other;
                }
            }
            1 => {
                if self.z(atm2) == 16
                    && self.z(atm3) == 16
                    && self.degree(atm2) == 2
                    && self.degree(atm3) == 2
                {
                    dl = compute14_dist_3d(bl1, bl2, bl3, ba12, ba23, PI / 2.0) - GEN_DIST_TOL;
                    du = dl + 2.0 * GEN_DIST_TOL;
                    ptype = Path14Type::Other;
                } else if self.check_macro_all_amide_ester_14(aid1, atm2, atm3, aid4)
                    || self.check_macro_all_amide_ester_14(aid4, atm3, atm2, aid1)
                {
                    // trans amide, plus the empirical +0.1
                    dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23) + 0.1;
                    ptype = Path14Type::Trans;
                    data.insert_trans(bid1, bid2, bid3);
                    du = dl;
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                } else if self.check_amide_ester_15(bid1, bid3, atm2, atm3)
                    || self.check_amide_ester_15(bid3, bid1, atm3, atm2)
                {
                    // FORCE_TRANS_AMIDES is undefined -> the #else branch: amide is cis, we're trans
                    if self.z(atm2) == 7
                        && self.degree(atm2) == 3
                        && self.z(aid1) == 1
                        && self.num_hs[atm2] == 1
                    {
                        set_the_bound = false;
                    } else {
                        dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        ptype = Path14Type::Trans;
                        data.insert_trans(bid1, bid2, bid3);
                    }
                    du = dl;
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                } else {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    ptype = Path14Type::Other;
                }
            }
            _ => {
                dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                ptype = Path14Type::Other;
            }
        }
        if set_the_bound {
            if (du - dl).abs() < DIST12_DELTA {
                dl -= GEN_DIST_TOL;
                du += GEN_DIST_TOL;
            }
            bounds.check_and_set(aid1, aid4, dl, du);
            data.paths14.push((bid1, bid2, bid3, ptype));
            let p = data.pid(aid1, aid4);
            data.visited14[p] = true;
        }
    }

    /// The `set14Bounds` driver: same-ring paths first, then every 1-4 bond path.
    fn run(&self, bounds: &mut BoundsMat, data: &mut ComputedData) {
        let nb = self.nb;
        let mut ring_bond_pairs: HashSet<(usize, usize)> = HashSet::new();
        let mut done_paths: HashSet<(usize, usize, usize)> = HashSet::new();
        let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();

        // first: 1-4 atoms in the same ring
        for ri in 0..self.bond_rings.len() {
            let bring = self.bond_rings[ri].clone();
            let rsize = bring.len();
            if rsize < 3 {
                continue;
            }
            let mut bid1 = bring[rsize - 1];
            for i in 0..rsize {
                let bid2 = bring[i];
                let bid3 = bring[(i + 1) % rsize];
                ring_bond_pairs.insert((bid1, bid2));
                ring_bond_pairs.insert((bid2, bid1));
                done_paths.insert((bid1, bid2, bid3));
                done_paths.insert((bid3, bid2, bid1));
                if rsize > 5 {
                    if self.macrocycle14 && rsize >= MIN_MACROCYCLE_RING_SIZE {
                        self.set_macrocycle_all_in_same_ring_14(bid1, bid2, bid3, bounds, data);
                        bid_is_macrocycle.insert(bid2);
                    } else {
                        self.set_in_ring_14(bid1, bid2, bid3, rsize as i64, bounds, data);
                    }
                } else {
                    self.record_14_path(bid1, bid2, bid3, data);
                }
                bid1 = bid2;
            }
        }

        // second: every 1-4 path through a central bond
        for bid2 in 0..nb {
            let (aid2, aid3) = self.mol.bonds[bid2];
            for &(bid1, _) in &self.ab[aid2] {
                if bid1 == bid2 {
                    continue;
                }
                for &(bid3, _) in &self.ab[aid3] {
                    if bid3 == bid2 {
                        continue;
                    }
                    if done_paths.contains(&(bid1, bid2, bid3))
                        || done_paths.contains(&(bid3, bid2, bid1))
                    {
                        continue;
                    }
                    let in_ring_pair = ring_bond_pairs.contains(&(bid1, bid2))
                        || ring_bond_pairs.contains(&(bid2, bid1))
                        || ring_bond_pairs.contains(&(bid2, bid3))
                        || ring_bond_pairs.contains(&(bid3, bid2));
                    if in_ring_pair {
                        if self.macrocycle14 && bid_is_macrocycle.contains(&bid2) {
                            self.set_macrocycle_two_in_same_ring_14(bid1, bid2, bid3, bounds, data);
                        } else {
                            self.set_two_in_same_ring_14(bid1, bid2, bid3, bounds, data);
                        }
                    } else if self.num_bond_rings[bid2] > 0
                        && (self.num_bond_rings[bid1] > 0 || self.num_bond_rings[bid3] > 0)
                    {
                        // two ring bonds in different rings -> treated like in-ring, ringSize 0
                        self.set_in_ring_14(bid1, bid2, bid3, 0, bounds, data);
                    } else if self.num_bond_rings[bid2] > 0 {
                        self.set_in_ring_14(bid1, bid2, bid3, 0, bounds, data);
                    } else {
                        self.set_chain_14(bid1, bid2, bid3, bounds, data);
                    }
                }
            }
        }
    }
}

/// 1-4 distance bounds. Reproduces `set14Bounds`.
pub fn set14(
    mol: &Perceived,
    bounds: &mut BoundsMat,
    data: &mut ComputedData,
    dmat: &[f64],
    macrocycle14: bool,
    force_trans_amides: bool,
) {
    let ctx = S14::new(mol, dmat, macrocycle14, force_trans_amides);
    ctx.run(bounds, data);
}

// ---------------------------------------------------------------------------------------------
// set15Bounds: 1-5 distance bounds, chained across two 1-4 paths sharing a central bond.
// ---------------------------------------------------------------------------------------------

/// Tolerance around a 1-5 distance (`DIST15_TOL`).
const DIST15_TOL: f64 = 0.08;

fn compute15_cis_cis(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() - d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 + a23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute13_dist(d14, d4, a34 - ang143)
}
fn compute15_cis_trans(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() - d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 + a23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute13_dist(d14, d4, a34 + ang143)
}
fn compute15_trans_trans(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() + d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 - a23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute13_dist(d14, d4, a34 + ang143)
}
fn compute15_trans_cis(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() + d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 - a23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute13_dist(d14, d4, a34 - ang143)
}

impl<'a> S14<'a> {
    /// `_set15BoundsHelper`: extend the 1-4 path `bid1-bid2-bid3` (of class `ptype`) by every bond
    /// off its far atom, setting the 1-5 bound from the two paths' cis/trans classes.
    #[allow(unused_assignments)]
    fn set15_helper(
        &self,
        bid1: usize,
        bid2: usize,
        bid3: usize,
        ptype: Path14Type,
        bounds: &mut BoundsMat,
        data: &mut ComputedData,
    ) {
        let na = self.n;
        let aid2 = data.bond_adj_atom(bid1, bid2);
        let aid1 = self.other_atom(bid1, aid2);
        let aid3 = data.bond_adj_atom(bid2, bid3);
        let aid4 = self.other_atom(bid3, aid3);
        let d1 = data.bond_lengths[bid1];
        let d2 = data.bond_lengths[bid2];
        let d3 = data.bond_lengths[bid3];
        let ang12 = data.bond_angle(bid1, bid2);
        let ang23 = data.bond_angle(bid2, bid3);

        for i in 0..self.nb {
            // bond i must extend the path off aid4 (share it with bid3)
            if data.bond_adj[bid3 * self.nb + i] != aid4 as i64 {
                continue;
            }
            let aid5 = self.other_atom(i, aid4);
            if data.visited_upto14(aid1, aid5) {
                return;
            }
            if self.dmat_dist(aid1, aid5) < 3.9 {
                continue;
            }
            if aid1 == aid5 {
                continue;
            }
            let p1 = aid1 * na + aid5;
            let p2 = aid5 * na + aid1;
            if bounds.lower(aid1, aid5) >= DIST12_DELTA
                && !data.set15_atoms[p1]
                && !data.set15_atoms[p2]
            {
                continue;
            }
            let d4 = data.bond_lengths[i];
            let ang34 = data.bond_angle(bid3, i);
            let path_id = data.path_key(bid2, bid3, i);
            let is_cis = data.cis_paths.contains(&path_id);
            let is_trans = data.trans_paths.contains(&path_id);
            let mut du = -1.0;
            let mut dl = 0.0;
            match ptype {
                Path14Type::Cis => {
                    if is_cis {
                        dl = compute15_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                        du = dl + DIST15_TOL;
                        dl -= DIST15_TOL;
                    } else if is_trans {
                        dl = compute15_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                        du = dl + DIST15_TOL;
                        dl -= DIST15_TOL;
                    } else {
                        dl = compute15_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                        du = compute15_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                    }
                }
                Path14Type::Trans => {
                    if is_cis {
                        dl = compute15_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                        du = dl + DIST15_TOL;
                        dl -= DIST15_TOL;
                    } else if is_trans {
                        dl = compute15_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                        du = dl + DIST15_TOL;
                        dl -= DIST15_TOL;
                    } else {
                        dl = compute15_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                        du =
                            compute15_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                    }
                }
                Path14Type::Other => {
                    if is_cis {
                        dl = compute15_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                        du = compute15_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                    } else if is_trans {
                        dl = compute15_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                        du =
                            compute15_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                    } else {
                        let vw1 = valence::rvdw(self.z(aid1)).unwrap_or(0.0);
                        let vw5 = valence::rvdw(self.z(aid5)).unwrap_or(0.0);
                        dl = VDW_SCALE_15 * (vw1 + vw5);
                    }
                }
            }
            if du < 0.0 {
                du = MAX_UPPER;
            }
            bounds.check_and_set(aid1, aid5, dl, du);
            data.set15_atoms[aid1 * na + aid5] = true;
            data.set15_atoms[aid5 * na + aid1] = true;
        }
    }
}

/// 1-5 distance bounds. Reproduces `set15Bounds`: for each recorded 1-4 path, extend it both ways.
pub fn set15(mol: &Perceived, bounds: &mut BoundsMat, data: &mut ComputedData, dmat: &[f64]) {
    let ctx = S14::new(mol, dmat, true, true);
    let paths = data.paths14.clone();
    for (bid1, bid2, bid3, ptype) in paths {
        ctx.set15_helper(bid1, bid2, bid3, ptype, bounds, data);
        ctx.set15_helper(bid3, bid2, bid1, ptype, bounds, data);
    }
}

/// van der Waals lower bounds for every pair not already constrained. Reproduces
/// `setLowerBoundVDW`: an H-bond donor-H / acceptor pair gets 1.8 Å; pairs 4 and 5 bonds apart get a
/// scaled sum of vdW radii; everything else the full sum. `dmat` is the topological distance matrix.
///
/// Faithful to a RDKit quirk: the outer loop runs `i` from 1, so atom 0's donor/acceptor flags are
/// never set and it is treated as neither — reproduced here rather than "corrected".
pub fn set_lower_bound_vdw(mol: &Perceived, bounds: &mut BoundsMat, dmat: &[f64]) {
    let n = mol.atomic_numbers.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b) in &mol.bonds {
        adj[a].push(b);
        adj[b].push(a);
    }
    let z = &mol.atomic_numbers;
    let is_acceptor = |i: usize| z[i] == 7 || z[i] == 8;
    let is_h_in_donor = |i: usize| z[i] == 1 && adj[i].iter().any(|&nb| z[nb] == 7 || z[nb] == 8);

    let mut hin = vec![false; n];
    let mut acc = vec![false; n];
    for i in 1..n {
        let vw1 = valence::rvdw(z[i]).unwrap_or(0.0);
        hin[i] = is_h_in_donor(i);
        acc[i] = is_acceptor(i);
        for j in 0..i {
            let vw2 = valence::rvdw(z[j]).unwrap_or(0.0);
            if bounds.lower(i, j) < DIST12_DELTA {
                let d = dmat[i * n + j];
                let lb = if (hin[i] && acc[j]) || (acc[i] && hin[j]) {
                    H_BOND_LENGTH
                } else if d == 4.0 {
                    VDW_SCALE_15 * (vw1 + vw2)
                } else if d == 5.0 {
                    (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2)
                } else {
                    vw1 + vw2
                };
                bounds.set_lower(i, j, lb);
            }
        }
    }
}

/// The topological bounds matrix with the stages selected by the flags — RDKit's `setTopolBounds`
/// with the matching `set13`/`set14`/`set15`. Returns the whole [`BoundsMat`] so callers can see its
/// `bad` flag (an infeasible matrix RDKit would have thrown on); the flat `data` matches RDKit's
/// `getData()`.
fn topol_bounds_mat(
    mol: &Perceived,
    set13_flag: bool,
    set14_flag: bool,
    set15_flag: bool,
) -> BoundsMat {
    let n = mol.atomic_numbers.len();
    let nb = mol.bonds.len();
    let mut bounds = BoundsMat::new(n);
    let mut data = ComputedData::new(n, nb);
    let dmat = dist::distance_matrix(n, &mol.bonds);
    set12(mol, &mut bounds, &mut data);
    if set13_flag {
        set13(mol, &mut bounds, &mut data);
    }
    if set14_flag {
        // setTopolBounds is called with useMacrocycle14config=true, forceTransAmides=true
        set14(mol, &mut bounds, &mut data, &dmat, true, true);
    }
    if set15_flag {
        set15(mol, &mut bounds, &mut data, &dmat);
    }
    set_lower_bound_vdw(mol, &mut bounds, &dmat);
    bounds
}

/// The topological bounds matrix with the stages selected by the flags; the flat result matches
/// RDKit's `getData()`. Ignores feasibility — for staged dumps and benchmarks that want the raw
/// matrix regardless; the spec path uses [`bounds_full_checked`].
pub fn topol_bounds(
    mol: &Perceived,
    set13_flag: bool,
    set14_flag: bool,
    set15_flag: bool,
) -> Vec<f64> {
    topol_bounds_mat(mol, set13_flag, set14_flag, set15_flag).data
}

/// The full `setTopolBounds` (all layers) — the raw, pre-smoothing bounds matrix.
pub fn bounds_full(mol: &Perceived) -> Vec<f64> {
    topol_bounds(mol, true, true, true)
}

/// The full `setTopolBounds`, rejecting an infeasible matrix the way RDKit's `_checkAndSetBounds`
/// invariants do (reported as "bad lower bound"). `Err` means the molecule is unbuildable.
pub fn bounds_full_checked(mol: &Perceived) -> Result<Vec<f64>, String> {
    let m = topol_bounds_mat(mol, true, true, true);
    if m.bad {
        return Err("bad lower bound".to_string());
    }
    Ok(m.data)
}

/// The bounds matrix with only 1-2 bounds and van der Waals lower bounds set.
pub fn bounds_set12_vdw(mol: &Perceived) -> Vec<f64> {
    topol_bounds(mol, false, false, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_layout() {
        let m = BoundsMat::new(3);
        let at = |i: usize, j: usize| m.data[i * 3 + j];
        assert_eq!(at(0, 1), MAX_UPPER);
        assert_eq!(at(1, 0), 0.0);
        assert_eq!(at(0, 0), 0.0);
    }

    #[test]
    fn upper_lower_addressing() {
        let mut m = BoundsMat::new(3);
        m.set_upper(2, 0, 5.0);
        m.set_lower(0, 2, 1.0);
        let at = |i: usize, j: usize| m.data[i * 3 + j];
        assert_eq!(at(0, 2), 5.0);
        assert_eq!(at(2, 0), 1.0);
        assert_eq!(m.lower(0, 2), 1.0);
    }

    #[test]
    fn compute13_is_law_of_cosines() {
        // a right angle: the 1-3 distance is the hypotenuse of legs d1, d2
        let d = compute13_dist(3.0, 4.0, PI / 2.0);
        assert!((d - 5.0).abs() < 1e-12);
    }
}
