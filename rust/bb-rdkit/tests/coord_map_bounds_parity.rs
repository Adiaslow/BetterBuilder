//! Parity gate for the core-pin bounds tightening (`bb_embed::bounds::coord_map_bounds`), the one
//! deterministic routine on the real build path (`embed_recipe` → `bb-output::assemble`) that was
//! previously only self-tested. RDKit tightens a coordMap'd bounds matrix in `setupInitialBoundsMatrix`:
//! raw topology bounds → `adjustBoundsMatFromCoordMap` (pin each pinned pair to its exact distance)
//! → `triangleSmoothBounds(tol=0.05)`. Native instead adjusts+re-smooths the already-smoothed
//! `spec.bounds`; this gate proves that shortcut yields the identical matrix. Pins are taken from a
//! real embedded conformer, so their mutual distances are geometrically realizable and smoothing
//! succeeds (RDKit's recompute-relaxed fallback — where native falls back to the un-pinned base — is
//! unreachable with conformer-derived pins, which is the only way the recipe ever produces pins).

/// Worst absolute element difference between native's f32 matrix and RDKit's f64 matrix.
fn max_abs_diff(native: &[f64], rdkit: &[f64]) -> f64 {
    assert_eq!(native.len(), rdkit.len(), "matrix size mismatch");
    native
        .iter()
        .zip(rdkit)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0f64, f64::max)
}

const MOLECULES: &[&str] = &[
    "OC(=O)c1ccccc1",                  // benzoic acid
    "C1CCC(CC1)C2CCCCC2",              // two rings
    "c1ccc2c(c1)cccc2C(=O)NCCc3ccccc3", // amide + aromatic
    "O=C1CCCCCCCCCCNC(=O)CCCCCCCCCN1",  // macrocycle — the real use case
    "C1CC2CCC1CC2",                    // bicyclooctane
];

#[test]
fn coord_map_bounds_match_rdkit() {
    let mut checked = 0usize;
    for &smiles in MOLECULES {
        let spec = bb_spec::build_native(smiles).expect("native spec");
        let n = spec.n_atoms;
        let conf = &bb_embed::embed(&spec, 1, 0xC0FFEE).expect("embed")[0].coords;

        // Pin roughly the first third of the atoms as a "core", at their embedded coordinates —
        // a geometrically consistent pin set (distances are realizable, so smoothing succeeds).
        let k = (n / 3).max(2).min(n);
        let mut pinned: Vec<Option<[f64; 3]>> = vec![None; n];
        let mut pin_idx: Vec<i32> = Vec::new();
        let mut pin_xyz: Vec<f64> = Vec::new();
        for a in 0..k {
            let p = [conf[a * 3], conf[a * 3 + 1], conf[a * 3 + 2]];
            pinned[a] = Some(p);
            pin_idx.push(a as i32);
            pin_xyz.extend_from_slice(&p);
        }

        let native = bb_embed::bounds::coord_map_bounds_f64(&spec.raw_bounds_f64, n, &pinned);
        let rdkit = bb_rdkit::coord_map_bounds(smiles, &pin_idx, &pin_xyz);

        assert!(
            !rdkit.is_empty(),
            "{smiles}: RDKit coordMap bounds empty (smoothing failed on conformer-derived pins?)"
        );
        assert_eq!(native.len(), n * n, "{smiles}: native matrix wrong size");

        let d = max_abs_diff(&native, &rdkit);
        // f64 throughout (raw_bounds_f64 + f64 coordMap) → machine precision, not f32.
        assert!(d < 1e-9, "{smiles}: coordMap bounds differ by {d:.2e}");
        checked += 1;
    }
    assert!(checked == MOLECULES.len(), "not all molecules checked");
}
