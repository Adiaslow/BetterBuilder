//! Parity gate for the embed accept/reject checks (`bb_embed::checks`) against RDKit's own
//! post-embed checks (`EmbeddingOps::checkTetrahedralCenters/checkChiralCenters/…`, plus the
//! planarity test inside `minimizeWithExpTorsions`). The `embed_checks` bridge oracle runs the
//! *real* RDKit routines on a supplied 3D conformer; here we feed native and RDKit the identical
//! conformer and require they make the identical accept/reject decision — not just on the embedded
//! (valid) geometry, but on reflected, jittered, and surgically-broken variants that push every
//! check across its threshold, so no check's agreement is vacuous.

use bb_embed::checks;

/// Deterministic LCG so the jitter is reproducible (no wall-clock / thread randomness).
struct Lcg(u64);
impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
    fn jitter(&mut self, c: &[f64], mag: f64) -> Vec<f64> {
        c.iter().map(|&x| x + mag * (self.next_f64() - 0.5)).collect()
    }
}

/// Native's six checks, in `passes_checks` order, as a bit vector (no pinned atoms — a free embed).
fn native_bits(spec: &bb_core::MoleculeSpec, c: &[f64]) -> Vec<i32> {
    let pinned = vec![None; spec.n_atoms];
    vec![
        checks::check_tetrahedral_centers(spec, c, &pinned) as i32,
        checks::check_chiral_centers(spec, c) as i32,
        checks::planarity_ok(spec, c) as i32,
        checks::double_bond_geometry_ok(spec, c) as i32,
        checks::final_chiral_checks(spec, c) as i32,
        checks::double_bond_stereo_ok(spec, c) as i32,
    ]
}

const NAMES: [&str; 6] = [
    "tetrahedral",
    "chiral",
    "planarity",
    "double-bond-geometry",
    "final-chiral",
    "double-bond-stereo",
];

/// Molecules chosen to exercise every check: chiral centers, potential-tetrahedral carbons
/// (bicyclooctane's ring carbons), cis/trans double bonds, and planar aromatic/amide systems.
const MOLECULES: &[&str] = &[
    "C[C@H](N)C(=O)O",        // alanine — one chiral center
    "O1C[C@]1(CCCCCCCC)[H]",  // epoxide + explicit-H chiral
    "C1CC[C@H]2CCCC[C@H]2C1", // cis-decalin — two ring chiral centers
    "C1CC2CCC1CC2",           // bicyclo[2.2.2]octane — potential-tetrahedral ring carbons
    "F/C=C/F",                // trans-difluoroethene — E double bond
    "F/C=C\\F",               // cis-difluoroethene — Z double bond
    "C/C=C/C=C/C",            // conjugated diene, both E
    "c1ccccc1C(=O)N/C=C/C",   // aromatic + amide + E enamide (planarity + stereo)
];

fn pt(c: &[f64], i: usize) -> [f64; 3] {
    [c[i * 3], c[i * 3 + 1], c[i * 3 + 2]]
}
fn set(c: &mut [f64], i: usize, p: [f64; 3]) {
    c[i * 3] = p[0];
    c[i * 3 + 1] = p[1];
    c[i * 3 + 2] = p[2];
}

#[test]
fn embed_checks_match_rdkit() {
    let mut rng = Lcg(0x51A7ED);
    let mut total_variants = 0usize;
    let mut fire = [0usize; 6];

    for &smiles in MOLECULES {
        let spec = bb_spec::build_native(smiles).expect("native spec");

        for conf in bb_embed::embed(&spec, 2, 0xC0FFEE) {
            let base = conf.coords;

            let mut variants: Vec<Vec<f64>> = vec![base.clone()];

            // Reflection: negate every atom's x → inverts every chiral volume (chiral/final fire).
            let mut reflected = base.clone();
            for x in reflected.iter_mut().step_by(3) {
                *x = -*x;
            }
            variants.push(reflected);

            // Jitters of growing magnitude walk planarity / final-chiral-bounds across their thresholds.
            for &mag in &[0.3f64, 0.7, 1.2, 2.0, 3.0] {
                variants.push(rng.jitter(&base, mag));
            }

            // Surgically flatten the first tetrahedral center: collapse its four substituents almost
            // onto the center → chiral volume ~0 → the tetrahedral volume test must fire.
            if let Some(ts) = spec.tetrahedral_centers.first() {
                let mut v = base.clone();
                let ctr = pt(&base, ts.center as usize);
                for &a in &ts.atoms {
                    let p = pt(&base, a as usize);
                    let near = [
                        ctr[0] + (p[0] - ctr[0]) * 0.02,
                        ctr[1] + (p[1] - ctr[1]) * 0.02,
                        ctr[2] + (p[2] - ctr[2]) * 0.02,
                    ];
                    set(&mut v, a as usize, near);
                }
                variants.push(v);
            }

            // Surgically linearize the first double-bond end (a0-a1-a2): put a0 = 2*a1 - a2 so the
            // three are collinear → doubleBondGeometryChecks' "not linear" test must fire.
            if let Some(&[e0, e1, e2]) = spec.double_bond_ends.first() {
                let mut v = base.clone();
                let (p1, p2) = (pt(&base, e1 as usize), pt(&base, e2 as usize));
                set(&mut v, e0 as usize, [2.0 * p1[0] - p2[0], 2.0 * p1[1] - p2[1], 2.0 * p1[2] - p2[2]]);
                variants.push(v);
            }

            for v in &variants {
                let rdkit = bb_rdkit::embed_checks(smiles, v);
                assert_eq!(rdkit.len(), 6, "{smiles}: oracle returned {} bits", rdkit.len());
                let native = native_bits(&spec, v);

                for k in 0..6 {
                    assert_eq!(
                        native[k], rdkit[k],
                        "{smiles}: check '{}' disagrees (native={}, rdkit={})",
                        NAMES[k], native[k], rdkit[k]
                    );
                    if rdkit[k] == 0 {
                        fire[k] += 1;
                    }
                }
                total_variants += 1;
            }
        }
    }

    // No check's agreement may be vacuous: every one of the six must have rejected at least once
    // (with native and RDKit agreeing on that rejection, per the assert above).
    for k in 0..6 {
        assert!(
            fire[k] > 0,
            "check '{}' never fired across {total_variants} variants — its parity is unproven",
            NAMES[k]
        );
    }
}
