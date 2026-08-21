//! bb-embed — the ETKDG macrocycle embed engine.
//!
//! Consumes a [`bb_core::MoleculeSpec`] and produces conformers by staged distance-geometry
//! minimization: init → Stage A (4D) → Stage B (4th-dimension squeeze) → Stage C (3D) → acceptance
//! checks with retry. Entry points: [`embed`] (independent conformers) and [`embed_recipe`] (the
//! two-stage core-pin recipe).

pub mod bounds;
pub mod checks;
pub mod forcefield;
pub mod init;
pub mod minimize;

use bb_core::MoleculeSpec;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;

/// `firstMinimization`'s per-atom energy reject threshold (RDKit `MAX_MINIMIZED_E_PER_ATOM`): a
/// Stage-A result with `calcEnergy()/nAtoms` at or above this is discarded and the conformer redrawn.
const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;

/// A generated conformer: `n_atoms * 3` xyz, row-major.
#[derive(Clone, Debug)]
pub struct Conformer {
    pub coords: Vec<f64>,
}

impl Conformer {
    #[inline]
    pub fn atom(&self, i: usize) -> [f64; 3] {
        [
            self.coords[3 * i],
            self.coords[3 * i + 1],
            self.coords[3 * i + 2],
        ]
    }
}

/// Initialization for the 4D Stage-A embed, matching RDKit's `useRandomCoords`.
#[derive(Clone, Copy, PartialEq)]
pub enum InitMode {
    /// Metric-matrix eigenvalue embedding (`useRandomCoords=False`) — RDKit's seed/core embeds.
    MetricMatrix,
    /// Random box (`useRandomCoords=True`) — RDKit's sidechain embeds.
    Random,
}

/// Generate one accepted conformer (3D, `n_atoms*3`): retry [`embed_attempt`] with fresh inits
/// until it passes [`checks::passes_checks`]. The attempt budget is `(10 × n_atoms).clamp(1, 24)`.
/// NOTE: RDKit's embedder uses `maxIterations = 10 * numAtoms` uncapped (`Embedder.cpp`), so for
/// n ≥ 3 native gives up sooner. Whether the 24 cap matches Divya's actual yield (the oracle — she
/// may set `maxIterations`) is a stochastic-layer question that can only be settled by the
/// end-to-end ensemble/yield comparison vs Divya, not per-unit vs RDKit; it is intentionally left at
/// the calibrated cap until that comparison is run. If none pass, the last attempt is returned.
///
/// `pinned` is empty (free embed) or length `n_atoms`; `Some(xyz)` holds that atom at `xyz`
/// throughout.
fn embed_one(
    spec: &MoleculeSpec,
    rng: &mut StdRng,
    pinned: &[Option<[f64; 3]>],
    init_mode: InitMode,
) -> Vec<f64> {
    let max_attempts = (10 * spec.n_atoms).clamp(1, 24);
    let mut last = Vec::new();
    for _ in 0..max_attempts {
        // `None` = a rejected attempt (init-level eigenvalue reject or the per-atom energy reject);
        // RDKit re-draws a fresh distance matrix on either, which the next loop iteration does.
        if let Some(c) = embed_attempt(spec, rng, pinned, init_mode) {
            if checks::passes_checks(spec, &c, pinned) {
                return c;
            }
            last = c;
        }
    }
    // Budget exhausted with nothing accepted: return the last completed (unaccepted) attempt, or an
    // empty conformer if every attempt was rejected (RDKit likewise yields no conformer here — the
    // ensemble consumers skip a conformer whose length != n*3). The 24-cap vs RDKit's uncapped 10·n
    // is a yield-calibration question for the Phase-3 ensemble comparison vs Divya.
    last
}

/// One embedding attempt (init → Stage A → Stage B → Stage C). Returns `None` on RDKit's embed-loop
/// rejects: the `computeInitialCoords` eigenvalue reject (via [`init::metric_matrix`]) or the
/// `firstMinimization` per-atom energy reject (`MAX_MINIMIZED_E_PER_ATOM = 0.05`) — either causes a
/// fresh re-draw. `Some(3D coords)` for a completed attempt (still subject to `passes_checks`).
fn embed_attempt(
    spec: &MoleculeSpec,
    rng: &mut StdRng,
    pinned: &[Option<[f64; 3]>],
    init_mode: InitMode,
) -> Option<Vec<f64>> {
    let n = spec.n_atoms;
    let fixed: Vec<bool> = if pinned.is_empty() {
        Vec::new()
    } else {
        pinned.iter().map(|p| p.is_some()).collect()
    };
    // RDKit sets basinThresh=1e8 (keep all distance terms) for random coords, 5.0 otherwise.
    let basin = match init_mode {
        InitMode::MetricMatrix => forcefield::BASIN_DEFAULT,
        InitMode::Random => forcefield::BASIN_ALL,
    };

    // Stage A init (4D), then pinned atoms placed at their frozen coords (4th dim 0). The metric-matrix
    // init can reject (degenerate distance matrix) — `?` propagates that as a re-draw.
    let mut a = match init_mode {
        InitMode::MetricMatrix => init::metric_matrix(spec, 4, rng)?,
        InitMode::Random => init::random_box(n, 4, 10.0, rng),
    };
    for (atom, p) in pinned.iter().enumerate() {
        if let Some(xyz) = p {
            a[atom * 4] = xyz[0];
            a[atom * 4 + 1] = xyz[1];
            a[atom * 4 + 2] = xyz[2];
            a[atom * 4 + 3] = 0.0;
        }
    }
    let a4 = minimize::minimize_stage_a(spec, a, 4, &fixed, basin, 400);

    // firstMinimization per-atom energy reject: RDKit sets gotCoords=false (→ re-draw) when
    // `calcEnergy()/nAtoms >= MAX_MINIMIZED_E_PER_ATOM (0.05)`. `stage_a_reject_energy` is RDKit's
    // calcEnergy convention (machine-precision-gated), so the decision matches exactly.
    if forcefield::stage_a_reject_energy(spec, &a4, 4, basin) >= MAX_MINIMIZED_E_PER_ATOM * n as f64 {
        return None;
    }

    // Stage B (4th-dim squeeze): RDKit runs `minimizeFourthDimension` when there are chiral centers
    // or we started from random coords (`chiralCenters>0 || useRandomCoords`).
    let b4 = if init_mode == InitMode::Random || !spec.chiral_sets.is_empty() {
        minimize::minimize_stage_b(spec, a4, 4, &fixed, basin, 200)
    } else {
        a4
    };

    // Project 4D -> 3D (pinned atoms keep their frozen xyz).
    let mut c3 = Vec::with_capacity(n * 3);
    for atom in 0..n {
        for c in 0..3 {
            c3.push(b4[atom * 4 + c]);
        }
    }

    // Stage C: constraints from the projected geometry, then torsions + impropers + constraints in 3D.
    let (dist_c, angle_c) = forcefield::build_stage_c_constraints(spec, &c3);
    Some(minimize::minimize_stage_c(spec, &dist_c, &angle_c, c3, &fixed, 300))
}

/// Embed `n_conf` independent conformers for `spec` from a single RNG seeded by `seed`. No
/// core-pinning; metric-matrix init.
pub fn embed(spec: &MoleculeSpec, n_conf: usize, seed: u64) -> Vec<Conformer> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n_conf)
        .map(|_| Conformer {
            coords: embed_one(spec, &mut rng, &[], InitMode::MetricMatrix),
        })
        .collect()
}

/// RNG seed for the conformer at core seed `j`, sidechain index `k`, derived from `base`.
#[inline]
fn mix_seed(base: u64, j: u64, k: u64) -> u64 {
    base ^ (j.wrapping_add(1)).wrapping_mul(0x9E3779B97F4A7C15)
        ^ (k.wrapping_add(1)).wrapping_mul(0xC2B2AE3D27D4EB4F)
}

/// The two-stage core-pin recipe. For each of `spec.core_seeds` seeds `j`: embed one free conformer
/// (metric-matrix init, RNG seed `base_seed + j`), hold `spec.pin_atoms` at those coordinates and
/// tighten the bounds from them ([`bounds::coord_map_bounds_f64`]), then embed `spec.sidechain_confs`
/// conformers (random init, one RNG seed per conformer) around the held core. Returns
/// `core_seeds × sidechain_confs` conformers, generated in parallel over rayon.
///
/// Falls back to 200 independent [`embed`] conformers when `core_seeds` or `sidechain_confs` is 0.
pub fn embed_recipe(spec: &MoleculeSpec, base_seed: u64) -> Vec<Conformer> {
    if spec.core_seeds == 0 || spec.sidechain_confs == 0 {
        return embed(spec, 200, base_seed);
    }
    let n = spec.n_atoms;
    let (n_seed, n_side) = (spec.core_seeds as usize, spec.sidechain_confs as usize);

    // Phase 1: core seed conformers — independent per j, embedded in parallel.
    let seed_coords: Vec<Vec<f64>> = (0..n_seed)
        .into_par_iter()
        .map(|j| {
            let mut rng = StdRng::seed_from_u64(base_seed + j as u64);
            embed_one(spec, &mut rng, &[], InitMode::MetricMatrix)
        })
        .collect();

    // Phase 2: per-seed setup — freeze the pin set + coordMap-tighten the bounds (serial, ~cheap).
    let per_seed: Vec<(Vec<Option<[f64; 3]>>, MoleculeSpec)> = seed_coords
        .iter()
        .map(|sc| {
            let mut pinned = vec![None; n];
            for &a in &spec.pin_atoms {
                let a = a as usize;
                pinned[a] = Some([sc[a * 3], sc[a * 3 + 1], sc[a * 3 + 2]]);
            }
            // Tighten the RAW bounds (f64, matching RDKit's coordMap path which starts pre-smoothing)
            // with the pinned distances; the sidechain embed reads bounds_f64 via ub64. Fall back to
            // the smoothed bounds for pre-field specs. `with_seed_bounds` shares the topology and skips
            // cloning the two raw-bounds matrices the sidechain embed never reads.
            let base = if spec.raw_bounds_f64.is_empty() { &spec.bounds_f64 } else { &spec.raw_bounds_f64 };
            let sc_spec = spec.with_seed_bounds(bounds::coord_map_bounds_f64(base, n, &pinned));
            (pinned, sc_spec)
        })
        .collect();

    // Phase 3: all sidechain conformers as one flat parallel task pool (best load balancing).
    let tasks: Vec<(usize, usize)> = (0..n_seed)
        .flat_map(|j| (0..n_side).map(move |k| (j, k)))
        .collect();
    tasks
        .into_par_iter()
        .map(|(j, k)| {
            let (pinned, sc_spec) = &per_seed[j];
            let mut rng = StdRng::seed_from_u64(mix_seed(base_seed, j as u64, k as u64));
            Conformer {
                coords: embed_one(sc_spec, &mut rng, pinned, InitMode::Random),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bb_core::MoleculeSpec;

    #[test]
    fn embed_produces_finite_conformers() {
        let n = 6;
        let mut bounds = vec![0.0f32; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                bounds[i * n + j] = 3.0;
                bounds[j * n + i] = 1.3;
            }
        }
        let spec = MoleculeSpec {
            n_atoms: n,
            dim: 4,
            bounds,
            ..Default::default()
        };
        let confs = embed(&spec, 3, 42);
        assert_eq!(confs.len(), 3);
        for c in &confs {
            assert_eq!(c.coords.len(), n * 3);
            assert!(c.coords.iter().all(|x| x.is_finite()));
        }
    }
}
