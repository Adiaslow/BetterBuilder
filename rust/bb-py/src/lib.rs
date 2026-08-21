//! Python bindings for BetterBuilder — a PyO3 wrapper over [`bb_spec::build_native`] (SMILES →
//! distance-geometry problem, RDKit-free) and [`bb_embed::embed_recipe`] (the ETKDG core-pin
//! recipe), exposed as one function:
//!
//! ```python
//! import betterbuilder
//! result = betterbuilder.embed("C1CCOCC1", seed=210185)
//! result.n_atoms          # heavy atoms + H, RDKit (SmilesToMol+AddHs) order
//! result.conformers[0]    # flat [x0, y0, z0, x1, …] for conformer 0
//! ```

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// A conformer ensemble for one molecule. Immutable; indexable/len over its conformers.
#[pyclass(frozen)]
struct Conformers {
    /// Atom count (heavy atoms + explicit H), in RDKit `SmilesToMol` + `AddHs` order.
    #[pyo3(get)]
    n_atoms: usize,
    /// One flat `n_atoms * 3` coordinate list per conformer.
    #[pyo3(get)]
    conformers: Vec<Vec<f64>>,
}

#[pymethods]
impl Conformers {
    fn __len__(&self) -> usize {
        self.conformers.len()
    }
    fn __repr__(&self) -> String {
        format!(
            "Conformers(n_atoms={}, n_conformers={})",
            self.n_atoms,
            self.conformers.len()
        )
    }
}

/// Embed a SMILES into a conformer ensemble via the ETKDG two-stage core-pin recipe. The conformer
/// count comes from the recipe counts in the spec, not from a caller argument.
///
/// Raises `ValueError` if the SMILES cannot be parsed.
#[pyfunction]
#[pyo3(signature = (smiles, seed = 210185))]
fn embed(smiles: &str, seed: u64) -> PyResult<Conformers> {
    let spec = bb_spec::build_native(smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let confs = bb_embed::embed_recipe(&spec, seed);
    Ok(Conformers {
        n_atoms: spec.n_atoms,
        conformers: confs.into_iter().map(|c| c.coords).collect(),
    })
}

/// BetterBuilder — macrocycle conformer generation.
#[pymodule]
fn betterbuilder(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(embed, m)?)?;
    m.add_class::<Conformers>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
