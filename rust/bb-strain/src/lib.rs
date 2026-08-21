//! bb-strain — UCSF torsion-strain energy, RDKit-free.
//!
//! A faithful port of the pipeline's `Torsion_Strain`/`TL_Functions`: match each rotatable torsion to
//! the most specific rule in the vendored Torsion Library ([`library`]), take the dihedral, and read
//! the strain penalty (an interpolated 36-bin histogram for `exact` rules, a `β₁·δ² + β₂·δ⁴` term for
//! `approximate` rules); the per-conformer **total** and **max** strain are what the db2 `S`-records
//! carry. SMARTS matching uses bb-perceive's engine (gated against RDKit), never RDKit at runtime.

pub mod compute;
pub mod library;

pub use compute::{match_library, strain};

#[cfg(test)]
mod tests {
    use super::library::{library, Method};

    #[test]
    fn library_parses() {
        let lib = library();
        // ~500 rules; a healthy mix of exact and approximate
        assert!(lib.len() > 400, "expected a full library, got {}", lib.len());
        let exact = lib.iter().filter(|r| matches!(r.method, Method::Exact { .. })).count();
        let approx = lib.len() - exact;
        assert!(exact > 0 && approx > 0, "exact {exact}, approx {approx}");
        // indices are strictly increasing in the lookup's specificity order
        assert!(lib.windows(2).all(|w| w[0].index < w[1].index));
        // every rule's SMARTS carries the four torsion atom maps
        assert!(lib.iter().all(|r| r.smarts.contains(":1") && r.smarts.contains(":4")));
    }
}
