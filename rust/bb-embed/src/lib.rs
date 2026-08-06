//! bb-embed — CPU reimplementation of the ETKDG macrocycle embed (the BetterBuilder candidate).
//!
//! Consumes a [`bb_core::MoleculeSpec`] (from `bb-rdkit`'s FFI to the patched RDKit) and produces
//! conformers by reimplementing ETKDG's staged distance-geometry minimization — the ~91% we
//! optimize. Numerics are best-in-class crates (argmin optimizer, rand); no Python. Validated by
//! ensemble RMSD vs the oracle (~2 Å bar), not by term parity.
//!
//! Status: first end-to-end path — random-box init → Stage-A minimize (distance + chiral + 4th-dim)
//! → 3D projection. Stages B/C (torsions/impropers/constraints) + the two-stage core-pin recipe are next.

pub mod bounds;
pub mod checks;
pub mod forcefield;
pub mod init;
pub mod minimize;

use bb_core::MoleculeSpec;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;

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

/// Generate one conformer (3D, `n_atoms*3`): 4D init (per `init_mode`) → **Stage A** (distance-
/// violation + chiral + 4th-dim) → project to 3D → **Stage C** (M6 torsions + impropers + distance
/// constraints).
///
/// `pinned` is empty (free embed) or length `n_atoms`; `Some(xyz)` freezes that atom to `xyz`
/// throughout (the candidate's coordMap/`fixedPoints` analog — used by the core-pin recipe).
/// Generate one accepted conformer: retry [`embed_attempt`] with fresh inits until it passes the
/// ETKDG acceptance checks ([`checks::passes_checks`]), matching RDKit's `embedPoints` retry loop.
/// Capped (RDKit uses `10×nAtoms`; we cap tighter for wall-clock — checks almost always pass on the
/// first try). If all attempts fail, the last is kept so the conformer count is preserved.
fn embed_one(
    spec: &MoleculeSpec,
    rng: &mut StdRng,
    pinned: &[Option<[f64; 3]>],
    init_mode: InitMode,
) -> Vec<f64> {
    let max_attempts = (10 * spec.n_atoms).clamp(1, 24);
    let mut last = Vec::new();
    for _ in 0..max_attempts {
        last = embed_attempt(spec, rng, pinned, init_mode);
        if checks::passes_checks(spec, &last, pinned) {
            return last;
        }
    }
    last
}

/// One embedding attempt (init → Stage A → Stage B → Stage C), no acceptance check.
fn embed_attempt(
    spec: &MoleculeSpec,
    rng: &mut StdRng,
    pinned: &[Option<[f64; 3]>],
    init_mode: InitMode,
) -> Vec<f64> {
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

    // Stage A init (4D), then pinned atoms placed at their frozen coords (4th dim 0).
    let mut a = match init_mode {
        InitMode::MetricMatrix => init::metric_matrix(spec, 4, rng),
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
    minimize::minimize_stage_c(spec, &dist_c, &angle_c, c3, &fixed, 300)
}

/// Embed `n_conf` independent conformers for `spec`, seeded by `seed` (no core-pinning). Uses the
/// metric-matrix init, matching RDKit `EmbedMolecule`'s `useRandomCoords=False` default.
pub fn embed(spec: &MoleculeSpec, n_conf: usize, seed: u64) -> Vec<Conformer> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n_conf)
        .map(|_| Conformer {
            coords: embed_one(spec, &mut rng, &[], InitMode::MetricMatrix),
        })
        .collect()
}

/// Deterministic per-conformer RNG seed, so every embed is independent (→ parallelizable) yet the
/// whole run stays reproducible for a given `base_seed`.
#[inline]
fn mix_seed(base: u64, j: u64, k: u64) -> u64 {
    base ^ (j.wrapping_add(1)).wrapping_mul(0x9E3779B97F4A7C15)
        ^ (k.wrapping_add(1)).wrapping_mul(0xC2B2AE3D27D4EB4F)
}

/// The faithful **two-stage core-pin recipe** (`build_ligands.py`), **parallelized with rayon** —
/// every conformer is independent, so the embed fans out across all cores. For each of `core_seeds`
/// seeds `j`, embed one free conformer (metric-matrix init, RNG `base_seed + j`), freeze the recipe's
/// `pin_atoms` to it + coordMap-tighten the bounds, then embed `sidechain_confs` conformers (random
/// init, per-conformer RNG) around the frozen core. Total = `core_seeds × sidechain_confs`. Falls back
/// to independent embeds if the spec carries no recipe data.
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
            let mut sc_spec = spec.clone();
            sc_spec.bounds = bounds::coord_map_bounds(&spec.bounds, n, &pinned);
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
