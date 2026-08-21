//! RDKit's per-element data, read from its own table.
//!
//! Implicit hydrogen counts follow the valence list rather than the octet rule or oxidation states:
//! `calcImplicitValence` takes the smallest entry that accommodates an atom's explicit valence, so
//! the hypervalent entries (S `[2, 4, 6]`, P `[3, 5]`, I `[1, 3, 5]`) decide the hydrogen count
//! wherever sulfones, phosphates or iodinanes appear. A valence of `-1` means unconstrained.
//!
//! The table is `rdkit-patch/atomic_data.tsv`, RDKit's `periodicTableAtomData` verbatim, parsed the
//! way RDKit parses it: fields are separated by runs of spaces and tabs with empty tokens dropped,
//! which matters because some rows leave the row-number column blank. Splitting on single tabs
//! instead shifts every later column. RDKit is BSD-3-Clause; see `rdkit-patch/README.md`.

use std::sync::OnceLock;

const ATOMIC_DATA: &str = include_str!("../../../rdkit-patch/atomic_data.tsv");

/// One row: the fields this crate uses.
#[derive(Debug, Clone, PartialEq)]
pub struct Element {
    pub atomic_number: u8,
    pub symbol: String,
    pub n_outer_elecs: u8,
    /// van der Waals radius (RDKit's `getRvdw`), used for VDW lower bounds and the missing-params
    /// bond-length fallback.
    pub rvdw: f64,
    /// Allowed valences, in the order RDKit lists them. `-1` means unconstrained.
    pub valences: Vec<i8>,
}

fn table() -> &'static Vec<Element> {
    static TABLE: OnceLock<Vec<Element>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut out = Vec::new();
        for line in ATOMIC_DATA.lines() {
            // split on runs of whitespace, dropping empties — boost::char_separator(" \t")
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 11 {
                continue;
            }
            let Ok(z) = f[0].parse::<u16>() else { continue };
            out.push(Element {
                atomic_number: u8::try_from(z).unwrap_or(0),
                symbol: f[1].to_string(),
                n_outer_elecs: f[7].parse().unwrap_or(0),
                // columns: anum symb row rCov rB0 rVdw mass nOuter isotope isotopeMass valences...
                rvdw: f[5].parse().unwrap_or(0.0),
                valences: f[10..].iter().filter_map(|v| v.parse().ok()).collect(),
            });
        }
        out.sort_by_key(|e| e.atomic_number);
        out
    })
}

/// The row for an atomic number, or `None` if the table has no such element.
pub fn element(z: u8) -> Option<&'static Element> {
    table().iter().find(|e| e.atomic_number == z && z > 0)
}

/// Allowed valences for an atomic number.
pub fn valence_list(z: u8) -> Option<&'static [i8]> {
    element(z).map(|e| e.valences.as_slice())
}

/// Outer-shell electron count, as RDKit's `getNouterElecs` reports it.
pub fn n_outer_elecs(z: u8) -> Option<u8> {
    element(z).map(|e| e.n_outer_elecs)
}

/// van der Waals radius, as RDKit's `getRvdw` reports it.
pub fn rvdw(z: u8) -> Option<f64> {
    element(z).map(|e| e.rvdw)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every value parsed here must match what RDKit's own accessors report. The fixture is
    /// captured by validation/parity/dump_valence.py; see its .provenance file.
    #[test]
    fn matches_rdkit_accessors() {
        let raw = include_str!("../../../validation/parity/fixtures/valence_rdkit.json");
        let ref_: serde_json::Value = serde_json::from_str(raw).expect("fixture parses");
        let mut checked = 0;
        for z in 1u8..=118 {
            let r = &ref_[z.to_string()];
            let e = element(z).unwrap_or_else(|| panic!("Z={z} missing from the table"));
            assert_eq!(e.symbol, r["symbol"].as_str().expect("symbol"), "Z={z}");
            let want: Vec<i8> = r["valence_list"]
                .as_array()
                .expect("valence_list")
                .iter()
                .map(|v| v.as_i64().expect("int") as i8)
                .collect();
            assert_eq!(e.valences, want, "valences for Z={z}");
            assert_eq!(
                i64::from(e.n_outer_elecs),
                r["n_outer_elecs"].as_i64().expect("outer"),
                "outer electrons for Z={z}"
            );
            assert!(
                (e.rvdw - r["rvdw"].as_f64().expect("rvdw")).abs() < 1e-12,
                "rvdw for Z={z}: {} vs {}",
                e.rvdw,
                r["rvdw"]
            );
            checked += 1;
        }
        assert_eq!(checked, 118);
    }

    #[test]
    fn symbols_agree_with_mendeleev() {
        for z in 1u8..=118 {
            let e = element(z).expect("element");
            assert_eq!(mendeleev::ALL_ELEMENTS[usize::from(z) - 1].symbol(), e.symbol);
        }
    }
}
