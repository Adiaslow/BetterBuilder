//! Construction of the distance-geometry bounds problem from a perceived molecule.
//!
//! This is the RDKit-free port of `setTopolBounds` (`GraphMol/DistGeomHelpers/BoundsMatrixBuilder`)
//! and the force-field parameters it reads. It consumes [`bb_perceive`]'s perception — elements,
//! rings, hybridization, valence — and produces the bounds a distance-geometry embedder needs.
//!
//! Kept separate from `bb-perceive` because these are force-field / geometry concerns, not molecular
//! perception (RDKit likewise keeps `ForceFieldHelpers` apart from `GraphMol` perception). The
//! dependency runs one way: bounds depend on perception.

pub mod bounds;
pub mod chirality;
pub mod dist;
pub mod uff;
