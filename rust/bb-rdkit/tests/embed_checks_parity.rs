//! Parity gate for the embed accept/reject checks (`bb_embed::checks`) against RDKit's own
//! post-embed checks (`EmbeddingOps::checkTetrahedralCenters/checkChiralCenters/…`, plus the
//! planarity test inside `minimizeWithExpTorsions`). The `embed_checks` bridge oracle runs the
//! *real* RDKit routines on a supplied 3D conformer; here we feed native and RDKit the identical
//! conformer and require they make the identical accept/reject decision — the behaviour the embed
//! acts on; which check rejects is not compared, since RDKit runs each check only after every earlier
//! one passed (`Embedder.cpp` `embedPoints`) and a rejected attempt is retried the same way whatever
//! rejected it. Variants go beyond the embedded (valid) geometry — reflected, jittered, and
//! surgically-broken — and the test requires an accepted variant (so a check that wrongly rejects
//! would be caught) and, for every check, a variant that check alone rejects in RDKit (so that check
//! wrongly accepting would flip the decision).

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

/// RDKit's six checks, in the order `embed_checks` reports them (its embedder's order).
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
    let mut accepted = 0usize;
    let mut decides = [0usize; 6];

    for &smiles in MOLECULES {
        let spec = bb_spec::build_native(smiles).expect("native spec");

        for conf in bb_embed::embed(&spec, 2, 0xC0FFEE).expect("embed") {
            let base = conf.coords;

            let mut variants: Vec<Vec<f64>> = vec![base.clone()];
            const REFLECTED: usize = 1; // index of the mirror image below

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

            // Flip the first stereo double bond b=c from E to Z (or Z to E) without distorting it: turn
            // the whole substituent on b's side 180 degrees about the b->c axis. Bond lengths and
            // angles are unchanged, so only the double-bond stereo check should reject. (Skipped when
            // b=c is in a ring: b's side then contains c.)
            if let Some(s) = spec.stereo_double_bonds.first() {
                let (b, c) = (s.atoms[1] as usize, s.atoms[2] as usize);
                let mut side = vec![false; spec.n_atoms];
                let mut stack = vec![b];
                side[b] = true;
                while let Some(a) = stack.pop() {
                    for &[i, j] in &spec.bonds {
                        let (i, j) = (i as usize, j as usize);
                        let next = if i == a { j } else if j == a { i } else { continue };
                        if side[next] || (a == b && next == c) {
                            continue;
                        }
                        side[next] = true;
                        stack.push(next);
                    }
                }
                if !side[c] {
                    let (pb, pc) = (pt(&base, b), pt(&base, c));
                    let axis = [pc[0] - pb[0], pc[1] - pb[1], pc[2] - pb[2]];
                    let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
                    let u = [axis[0] / len, axis[1] / len, axis[2] / len];
                    let mut v = base.clone();
                    for a in (0..spec.n_atoms).filter(|&a| side[a] && a != b) {
                        // 180-degree turn about the axis = reflection through the axis line
                        let p = pt(&base, a);
                        let t = (p[0] - pb[0]) * u[0] + (p[1] - pb[1]) * u[1] + (p[2] - pb[2]) * u[2];
                        let foot = [pb[0] + t * u[0], pb[1] + t * u[1], pb[2] + t * u[2]];
                        set(&mut v, a, [2.0 * foot[0] - p[0], 2.0 * foot[1] - p[1], 2.0 * foot[2] - p[2]]);
                    }
                    variants.push(v);
                }
            }

            for (vi, v) in variants.iter().enumerate() {
                let rdkit = bb_rdkit::embed_checks(smiles, v);
                assert_eq!(rdkit.len(), 6, "{smiles}: oracle returned {} bits", rdkit.len());
                // RDKit accepts the conformer only if every check passes; native's decision is the
                // one the embed uses (no pinned atoms — a free embed).
                let rdkit_accepts = rdkit.iter().all(|&b| b == 1);
                let native_accepts = checks::passes_checks(&spec, v, &[]);
                assert_eq!(
                    native_accepts, rdkit_accepts,
                    "{smiles} variant {vi}: native accepts {native_accepts}, RDKit accepts {rdkit_accepts} (RDKit checks {rdkit:?})"
                );
                accepted += usize::from(rdkit_accepts);
                match (0..6).filter(|&k| rdkit[k] == 0).collect::<Vec<_>>()[..] {
                    [k] => decides[k] += 1,
                    // RDKit's final chiral check starts by re-running the chiral check
                    // (`Embedder.cpp` `finalChiralChecks`), so a chiral failure always fails both.
                    // On the mirror image every distance and every centre-in-volume relation is
                    // preserved, so the final check's own tests pass and the chiral check alone decides.
                    [1, 4] if vi == REFLECTED => decides[1] += 1,
                    _ => {}
                }
                total_variants += 1;
            }
        }
    }

    // No agreement may be vacuous: some variant must be accepted, and every check must alone decide
    // some variant — there the decision turns on that check.
    assert!(accepted > 0, "no variant of {total_variants} was accepted — a wrongly rejecting check would go unseen");
    for (k, &n) in decides.iter().enumerate() {
        assert!(n > 0, "check '{}' never alone decided any of {total_variants} variants — its parity is unproven", NAMES[k]);
    }
}
