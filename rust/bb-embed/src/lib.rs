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

/// The safeguard against an input that cannot be embedded: [`embed`] makes at most this many times
/// as many draws as conformers requested, and [`embed_recipe`] uses at most this many times as many
/// cores as blocks requested, before reporting a [`Shortfall`].
pub const MAX_TRIES_PER_REQUEST: usize = 2;

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

/// One block of the core-pin recipe: the sidechain conformers embedded around one core, every one
/// holding `spec.pin_atoms` at that core's coordinates.
#[derive(Clone, Debug)]
pub struct Block {
    /// The core the block is built on, as numbered by [`recipe_core`].
    pub core: usize,
    pub conformers: Vec<Conformer>,
}

/// Fewer conformers than requested could be embedded, even after every redraw or replacement core
/// the entry point allows.
#[derive(Debug, thiserror::Error)]
#[error("embedded {embedded} of {requested} requested conformers")]
pub struct Shortfall {
    pub requested: usize,
    pub embedded: usize,
}

/// Why [`embed_recipe`] returned no ensemble. The first two mean the spec lacks an input the recipe
/// needs, as a spec written before those fields existed does, or the RDKit bridge spec, which is built
/// for comparison and never embedded; `bb_spec::build_native` supplies both.
#[derive(Debug, thiserror::Error)]
pub enum RecipeError {
    /// `core_seeds` or `sidechain_confs` is 0, so there is no recipe to run.
    #[error("the spec has no core-pin recipe (core_seeds {core_seeds}, sidechain_confs {sidechain_confs})")]
    NoRecipeCounts { core_seeds: u32, sidechain_confs: u32 },
    /// `raw_bounds_f64` is empty: the pre-smoothing bounds that RDKit's coordMap path tightens to a
    /// core's pins.
    #[error("the spec has no raw bounds matrix for the core-pin recipe to tighten")]
    NoRawBounds,
    #[error(transparent)]
    Shortfall(#[from] Shortfall),
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
/// until it passes [`checks::passes_checks`]. The attempt budget is RDKit's: `maxIterations = 0`,
/// which `build_ligands.py` leaves at its default, means `10 × numAtoms` (`Embedder.cpp`,
/// `embedPoints`). `None` if no attempt passes, as RDKit then yields no conformer.
///
/// `pinned` is empty (free embed) or length `n_atoms`; `Some(xyz)` holds that atom at `xyz`
/// throughout.
fn embed_one(
    spec: &MoleculeSpec,
    rng: &mut StdRng,
    pinned: &[Option<[f64; 3]>],
    init_mode: InitMode,
) -> Option<Vec<f64>> {
    let max_attempts = 10 * spec.n_atoms;
    for _ in 0..max_attempts {
        // `None` = a rejected attempt (init-level eigenvalue reject or the per-atom energy reject);
        // RDKit re-draws a fresh distance matrix on either, which the next loop iteration does.
        if let Some(c) = embed_attempt(spec, rng, pinned, init_mode) {
            if checks::passes_checks(spec, &c, pinned) {
                return Some(c);
            }
        }
    }
    None
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
/// core-pinning; metric-matrix init. A conformer that cannot be embedded is drawn again from the
/// same RNG, up to [`MAX_TRIES_PER_REQUEST`]` × n_conf` draws in all; past that the input is taken as
/// unembeddable and the result is a [`Shortfall`].
pub fn embed(spec: &MoleculeSpec, n_conf: usize, seed: u64) -> Result<Vec<Conformer>, Shortfall> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut conformers = Vec::with_capacity(n_conf);
    let mut draws_left = MAX_TRIES_PER_REQUEST * n_conf;
    while conformers.len() < n_conf && draws_left > 0 {
        draws_left -= 1;
        if let Some(coords) = embed_one(spec, &mut rng, &[], InitMode::MetricMatrix) {
            conformers.push(Conformer { coords });
        }
    }
    if conformers.len() == n_conf {
        Ok(conformers)
    } else {
        Err(Shortfall { requested: n_conf, embedded: conformers.len() })
    }
}

/// RNG seed of the recipe's core `j`, derived from `base`.
#[inline]
fn core_seed(base: u64, j: usize) -> u64 {
    base.wrapping_add(j as u64)
}

/// RNG seed of sidechain conformer `k` around core `j`, derived from `base`.
#[inline]
fn mix_seed(base: u64, j: usize, k: usize) -> u64 {
    base ^ (j as u64).wrapping_add(1).wrapping_mul(0x9E3779B97F4A7C15)
        ^ (k as u64).wrapping_add(1).wrapping_mul(0xC2B2AE3D27D4EB4F)
}

/// The recipe's core `j` for `base_seed`: one free conformer (metric-matrix init, nothing pinned), or
/// `None` if it cannot be embedded. [`embed_recipe`] builds block `s` around core `s`, and a
/// replacement block around the next unused core.
pub fn recipe_core(spec: &MoleculeSpec, base_seed: u64, j: usize) -> Option<Conformer> {
    let mut rng = StdRng::seed_from_u64(core_seed(base_seed, j));
    embed_one(spec, &mut rng, &[], InitMode::MetricMatrix).map(|coords| Conformer { coords })
}

/// [`embed_recipe`]'s block on core `j`, its `n_side` sidechains embedded in parallel. `None` if the
/// core or any of its sidechains cannot be embedded; the block's remaining sidechains are then skipped.
fn core_block(spec: &MoleculeSpec, base_seed: u64, j: usize, n_side: usize) -> Option<Block> {
    let n = spec.n_atoms;
    let core = recipe_core(spec, base_seed, j)?;
    let mut pinned = vec![None; n];
    for &a in &spec.pin_atoms {
        pinned[a as usize] = Some(core.atom(a as usize));
    }
    // Tighten the RAW bounds (f64, matching RDKit's coordMap path which starts pre-smoothing) with
    // the pinned distances; the sidechain embed reads bounds_f64 via ub64. `with_seed_bounds` shares
    // the topology and skips cloning the two raw-bounds matrices the sidechain embed never reads.
    let sc_spec = spec.with_seed_bounds(bounds::coord_map_bounds_f64(&spec.raw_bounds_f64, n, &pinned));
    let conformers = (0..n_side)
        .into_par_iter()
        .map(|k| {
            let mut rng = StdRng::seed_from_u64(mix_seed(base_seed, j, k));
            embed_one(&sc_spec, &mut rng, &pinned, InitMode::Random).map(|coords| Conformer { coords })
        })
        .collect::<Option<Vec<_>>>()?;
    Some(Block { core: j, conformers })
}

/// The two-stage core-pin recipe: `spec.core_seeds` blocks of `spec.sidechain_confs` conformers, in
/// block order, generated in parallel over rayon. Block `s` is built on core `s` ([`recipe_core`]):
/// `spec.pin_atoms` are held at that core's coordinates, the raw bounds are tightened to them
/// ([`bounds::coord_map_bounds_f64`]), and the sidechain conformers are embedded around the held core
/// (random init, one RNG seed per conformer).
///
/// Returns exactly `core_seeds` complete blocks, every conformer accepted. A block that cannot be
/// completed (its core, or a sidechain around it, fails to embed) is rebuilt around the next unused
/// core, `core_seeds`, `core_seeds + 1`, …, up to [`MAX_TRIES_PER_REQUEST`]` × core_seeds` cores in
/// all; past that the input is taken as unembeddable and the result is a [`Shortfall`]. This departs
/// from `build_ligands.py`, which leaves a failed core's conformers out.
pub fn embed_recipe(spec: &MoleculeSpec, base_seed: u64) -> Result<Vec<Block>, RecipeError> {
    if spec.core_seeds == 0 || spec.sidechain_confs == 0 {
        return Err(RecipeError::NoRecipeCounts {
            core_seeds: spec.core_seeds,
            sidechain_confs: spec.sidechain_confs,
        });
    }
    if spec.raw_bounds_f64.is_empty() {
        return Err(RecipeError::NoRawBounds);
    }
    let (n_seed, n_side) = (spec.core_seeds as usize, spec.sidechain_confs as usize);

    let mut blocks: Vec<Option<Block>> = vec![None; n_seed];
    let mut core_of: Vec<usize> = (0..n_seed).collect();
    let mut next_core = n_seed;
    let mut pending: Vec<usize> = (0..n_seed).collect();
    while !pending.is_empty() {
        let built: Vec<Option<Block>> = pending
            .par_iter()
            .map(|&s| core_block(spec, base_seed, core_of[s], n_side))
            .collect();
        let mut retry = Vec::new();
        for (s, block) in pending.into_iter().zip(built) {
            if block.is_some() {
                blocks[s] = block;
            } else if next_core < MAX_TRIES_PER_REQUEST * n_seed {
                core_of[s] = next_core;
                next_core += 1;
                retry.push(s);
            }
        }
        pending = retry;
    }

    let built: Vec<Block> = blocks.into_iter().flatten().collect();
    if built.len() == n_seed {
        Ok(built)
    } else {
        Err(Shortfall { requested: n_seed * n_side, embedded: built.len() * n_side }.into())
    }
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
        let confs = embed(&spec, 3, 42).expect("a realisable 6-atom spec embeds");
        assert_eq!(confs.len(), 3);
        for c in &confs {
            assert_eq!(c.coords.len(), n * 3);
            assert!(c.coords.iter().all(|x| x.is_finite()));
        }
    }
}
