//! Staged minimization: argmin L-BFGS (MoreThuente line search) over the Stage-A/B/C objectives
//! from [`crate::forcefield`].

use std::cell::RefCell;
use std::sync::OnceLock;

use argmin::core::{CostFunction, Error, Executor, Gradient};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use bb_core::MoleculeSpec;

use crate::forcefield::{
    build_dist_pairs, stage_a_energy_grad_pairs, stage_b_energy_grad_pairs, stage_c_energy_grad,
    AngleConstraint, DistConstraint, DistPair,
};

/// L-BFGS convergence tolerance on the raw gradient L2 norm (argmin `tolerance_grad`). Read once
/// from `BB_GRAD_TOL` if set, otherwise this default.
const GRAD_TOL_DEFAULT: f64 = 1e-2;
fn grad_tol() -> f64 {
    static T: OnceLock<f64> = OnceLock::new();
    *T.get_or_init(|| {
        std::env::var("BB_GRAD_TOL")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(GRAD_TOL_DEFAULT)
    })
}

/// Zero the gradient of fixed atoms so the optimizer does not move them. `fixed` is empty (no fixed
/// atoms) or length `n_atoms`.
#[inline]
fn zero_fixed(g: &mut [f64], fixed: &[bool], dim: usize) {
    for (a, &f) in fixed.iter().enumerate() {
        if f {
            for c in 0..dim {
                g[a * dim + c] = 0.0;
            }
        }
    }
}

/// Which staged objective a [`Problem`] evaluates.
enum Stage<'a> {
    A {
        basin: f64,
    },
    B {
        basin: f64,
    },
    C {
        dist_c: &'a [DistConstraint],
        angle_c: &'a [AngleConstraint],
    },
}

/// Memoized `(evaluated point, energy, gradient)` for one force-field evaluation.
type EvalCache = Option<(Vec<f64>, f64, Vec<f64>)>;

/// One minimization problem for argmin. Computes energy and gradient together and memoizes the last
/// evaluated point, so `cost(x)` and `gradient(x)` at the same `x` evaluate the field once. For
/// Stage A/B the basin-surviving distance pairs are pre-built once here (constant across the whole
/// minimization) instead of re-scanning all n²/2 pairs per eval.
struct Problem<'a> {
    spec: &'a MoleculeSpec,
    dim: usize,
    fixed: &'a [bool],
    stage: Stage<'a>,
    dist_pairs: Vec<DistPair>, // basin-surviving pairs for Stage A/B; empty for Stage C
    cache: RefCell<EvalCache>,
}

impl<'a> Problem<'a> {
    fn new(spec: &'a MoleculeSpec, dim: usize, fixed: &'a [bool], stage: Stage<'a>) -> Self {
        let dist_pairs = match &stage {
            Stage::A { basin } | Stage::B { basin } => build_dist_pairs(spec, *basin),
            Stage::C { .. } => Vec::new(),
        };
        Problem {
            spec,
            dim,
            fixed,
            stage,
            dist_pairs,
            cache: RefCell::new(None),
        }
    }

    fn compute(&self, p: &[f64]) -> (f64, Vec<f64>) {
        let (e, mut g) = match &self.stage {
            Stage::A { .. } => stage_a_energy_grad_pairs(self.spec, p, self.dim, &self.dist_pairs),
            Stage::B { .. } => stage_b_energy_grad_pairs(self.spec, p, self.dim, &self.dist_pairs),
            Stage::C { dist_c, angle_c } => stage_c_energy_grad(self.spec, dist_c, angle_c, p),
        };
        zero_fixed(&mut g, self.fixed, self.dim);
        (e, g)
    }

    /// Energy + gradient at `p`, reusing the cached result when `p` is the last point evaluated.
    fn eval(&self, p: &Vec<f64>) -> (f64, Vec<f64>) {
        if let Some((cp, ce, cg)) = self.cache.borrow().as_ref() {
            if cp == p {
                return (*ce, cg.clone());
            }
        }
        let (e, g) = self.compute(p);
        *self.cache.borrow_mut() = Some((p.clone(), e, g.clone()));
        (e, g)
    }
}

impl CostFunction for Problem<'_> {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        Ok(self.eval(p).0)
    }
}
impl Gradient for Problem<'_> {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, Error> {
        Ok(self.eval(p).1)
    }
}

/// Run L-BFGS (MoreThuente line search) on a staged [`Problem`] from `init`; returns the minimized
/// coordinates, or `init` unchanged if the optimizer returns an error.
fn run_min(problem: Problem, init: Vec<f64>, max_iters: u64) -> Vec<f64> {
    let solver = LBFGS::new(MoreThuenteLineSearch::new(), 7)
        .with_tolerance_grad(grad_tol())
        .expect("tolerance_grad >= 0");
    match Executor::new(problem, solver)
        .configure(|s| s.param(init.clone()).max_iters(max_iters))
        .run()
    {
        Ok(r) => r.state().best_param.clone().unwrap_or(init),
        Err(_) => init, // optimizer failure (rare; e.g. line-search stall) → keep the init
    }
}

/// Minimize the Stage-A objective from `init` (`n_atoms*dim` coords). `fixed` (empty or length
/// `n_atoms`) marks atoms held in place.
pub fn minimize_stage_a(
    spec: &MoleculeSpec,
    init: Vec<f64>,
    dim: usize,
    fixed: &[bool],
    basin: f64,
    max_iters: u64,
) -> Vec<f64> {
    run_min(
        Problem::new(spec, dim, fixed, Stage::A { basin }),
        init,
        max_iters,
    )
}

/// Minimize the Stage-B objective (4D 4th-dim squeeze) from `init`. `fixed` marks atoms held in place.
pub fn minimize_stage_b(
    spec: &MoleculeSpec,
    init: Vec<f64>,
    dim: usize,
    fixed: &[bool],
    basin: f64,
    max_iters: u64,
) -> Vec<f64> {
    run_min(
        Problem::new(spec, dim, fixed, Stage::B { basin }),
        init,
        max_iters,
    )
}

/// Minimize the Stage-C objective (3D) from `init`, using constraints pre-built from `init`'s
/// geometry (see [`crate::forcefield::build_stage_c_constraints`]). `fixed` marks atoms held in place.
pub fn minimize_stage_c(
    spec: &MoleculeSpec,
    dist_c: &[DistConstraint],
    angle_c: &[AngleConstraint],
    init: Vec<f64>,
    fixed: &[bool],
    max_iters: u64,
) -> Vec<f64> {
    run_min(
        Problem::new(spec, 3, fixed, Stage::C { dist_c, angle_c }),
        init,
        max_iters,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::{stage_a_energy_grad, BASIN_ALL};
    use crate::init::random_box;
    use bb_core::{ChiralSet, MoleculeSpec};
    use rand::SeedableRng;

    fn synthetic_spec() -> MoleculeSpec {
        let n = 6;
        let mut bounds = vec![0.0f32; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                bounds[i * n + j] = 3.0;
                bounds[j * n + i] = 1.3;
            }
        }
        MoleculeSpec {
            n_atoms: n,
            dim: 4,
            bounds,
            chiral_sets: vec![ChiralSet {
                center: 0,
                atoms: [1, 2, 3, 4],
                vol_lo: 5.0,
                vol_hi: 100.0,
                fused_small_rings: false,
            }],
            ..Default::default()
        }
    }

    #[test]
    fn minimize_lowers_energy() {
        let spec = synthetic_spec();
        let dim = 4;
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let init = random_box(spec.n_atoms, dim, 10.0, &mut rng);
        let e0 = stage_a_energy_grad(&spec, &init, dim, BASIN_ALL).0;
        let mini = minimize_stage_a(&spec, init, dim, &[], BASIN_ALL, 400);
        let e1 = stage_a_energy_grad(&spec, &mini, dim, BASIN_ALL).0;
        assert!(e1 < e0, "energy did not decrease: {e0:.3} -> {e1:.3}");
        // a satisfiable random spec should reach near-zero Stage-A energy per atom
        assert!(
            (e1 / spec.n_atoms as f64) < 0.05,
            "E/nAtoms too high: {:.4}",
            e1 / spec.n_atoms as f64
        );
    }
}
