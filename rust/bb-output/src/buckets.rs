//! Position clustering — a faithful port of `buckets2.py`. For one atom's position in each input
//! conformation, group the conformations whose positions coincide within `tolerance` into clusters,
//! using a spatial hash (bucket edge `tolerance/√3`, the largest cube inscribed in a tolerance-radius
//! sphere) plus a Moore neighborhood and position-dependent extra faces. Each cluster is a set of
//! conf indices; `bucket` accumulates them into `conf_clusters` keyed by that (sorted) set, in the
//! Python dict's insertion order — which downstream conf/set numbering depends on.

use std::collections::HashMap;

fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    let (dx, dy, dz) = (b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    dx * dx + dy * dy + dz * dz
}

/// Value stored per distinct position-sharing set: (assigned conf number, atom ids, their xyz).
pub type ClusterVal = (usize, Vec<usize>, Vec<[f64; 3]>);

pub struct Buckets {
    tolerance2: f64,
    bucketsize: f64,
    extrawidth: f64,
    moore: Vec<(i64, i64, i64)>,
    faces: [Vec<(i64, i64, i64)>; 6],
}

impl Buckets {
    pub fn new(tolerance: f64) -> Self {
        let bucketsize = tolerance / 3.0_f64.sqrt();
        let mut moore = Vec::new();
        for x in -1..=1 {
            for y in -1..=1 {
                for z in -1..=1 {
                    if (x, y, z) != (0, 0, 0) {
                        moore.push((x, y, z));
                    }
                }
            }
        }
        let face = |sel: fn(i64, i64) -> (i64, i64, i64)| {
            let mut v = Vec::new();
            for i in -1..=1 {
                for j in -1..=1 {
                    v.push(sel(i, j));
                }
            }
            v
        };
        let faces = [
            face(|i, j| (-2, i, j)),
            face(|i, j| (2, i, j)),
            face(|i, j| (i, -2, j)),
            face(|i, j| (i, 2, j)),
            face(|i, j| (i, j, -2)),
            face(|i, j| (i, j, 2)),
        ];
        Buckets { tolerance2: tolerance * tolerance, bucketsize, extrawidth: tolerance - bucketsize, moore, faces }
    }

    fn moore_faces(&self, xyz: [f64; 3], b: (i64, i64, i64)) -> Vec<usize> {
        let xd = xyz[0] - b.0 as f64 * self.bucketsize;
        let yd = xyz[1] - b.1 as f64 * self.bucketsize;
        let zd = xyz[2] - b.2 as f64 * self.bucketsize;
        let mut f = Vec::new();
        if xd < self.extrawidth {
            f.push(0);
        }
        if xd > self.bucketsize - self.extrawidth {
            f.push(1);
        }
        if yd < self.extrawidth {
            f.push(2);
        }
        if yd > self.bucketsize - self.extrawidth {
            f.push(3);
        }
        if zd < self.extrawidth {
            f.push(4);
        }
        if zd > self.bucketsize - self.extrawidth {
            f.push(5);
        }
        f
    }

    /// Cluster one atom's positions across confs. `xyz_data[conf]` is this atom's coordinate in each
    /// input conformation. Updates `conf_clusters` (insertion-ordered via `order`) and returns the
    /// number of distinct positions and the advanced `conf_num`.
    #[allow(clippy::too_many_arguments)]
    pub fn bucket(
        &self,
        xyz_data: &[[f64; 3]],
        atom_id: usize,
        mut conf_num: usize,
        clusters: &mut HashMap<Vec<usize>, ClusterVal>,
        order: &mut Vec<Vec<usize>>,
    ) -> (usize, usize) {
        // spatial hash: bucket -> (visited, conf indices), insertion-ordered
        let mut bucket_order: Vec<(i64, i64, i64)> = Vec::new();
        let mut buckets: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        let mut visited: HashMap<(i64, i64, i64), bool> = HashMap::new();
        for (conf, xyz) in xyz_data.iter().enumerate() {
            let h = (
                (xyz[0] / self.bucketsize).floor() as i64,
                (xyz[1] / self.bucketsize).floor() as i64,
                (xyz[2] / self.bucketsize).floor() as i64,
            );
            if let std::collections::hash_map::Entry::Vacant(e) = buckets.entry(h) {
                e.insert(vec![conf]);
                visited.insert(h, false);
                bucket_order.push(h);
            } else {
                buckets.get_mut(&h).unwrap().push(conf);
            }
        }

        let mut npos = 0;
        for h in &bucket_order {
            // each bucket seeds a cluster once; skip only if a prior seed absorbed all its confs
            if buckets[h].is_empty() {
                continue;
            }
            let seed_conf = buckets[h][0];
            let xyz = xyz_data[seed_conf];
            let mut curr: Vec<usize> = buckets[h].clone();

            let faces = self.moore_faces(xyz, *h);
            let mut locs: Vec<(i64, i64, i64)> = self.moore.iter().map(|n| (h.0 + n.0, h.1 + n.1, h.2 + n.2)).collect();
            for &fi in &faces {
                for n in &self.faces[fi] {
                    locs.push((h.0 + n.0, h.1 + n.1, h.2 + n.2));
                }
            }
            for loc in locs {
                // absorb a neighbour bucket's within-tolerance confs (skip missing / already-seeded)
                if !buckets.contains_key(&loc) || visited[&loc] {
                    continue;
                }
                let neigh = std::mem::take(buckets.get_mut(&loc).unwrap());
                let mut kept = Vec::new();
                for ni in neigh {
                    if dist2(xyz, xyz_data[ni]) <= self.tolerance2 {
                        curr.push(ni);
                    } else {
                        kept.push(ni);
                    }
                }
                *buckets.get_mut(&loc).unwrap() = kept;
            }
            visited.insert(*h, true);

            curr.sort_unstable();
            match clusters.get_mut(&curr) {
                Some((_, atomlist, xyzlist)) => {
                    atomlist.push(atom_id);
                    xyzlist.push(xyz);
                }
                None => {
                    clusters.insert(curr.clone(), (conf_num, vec![atom_id], vec![xyz]));
                    order.push(curr);
                    conf_num += 1;
                }
            }
            npos += 1;
        }
        (npos, conf_num)
    }
}
