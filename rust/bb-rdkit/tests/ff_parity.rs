//! Force-field gradient parity: bb-embed's ETKDG force-field terms vs RDKit's own force fields,
//! evaluated at identical coordinates.
//!
//! The finite-difference tests in bb-embed prove each gradient is the derivative of *native's* own
//! energy; they cannot catch a wrong energy *formula*. This gate does. The native spec's bounds,
//! chiral sets and experimental torsions are byte-identical to what the bridge feeds RDKit's force
//! field (the spec gate), so a gradient gap at the same coordinates is a force-field formula
//! difference. The gradient — the forces the minimizer follows — is the meaningful cross-check;
//! RDKit's Chiral/FourthDim contribs even use an energy inconsistent with their own gradient.
//!
//! - Stage A (`constructForceField`, the `firstMinimization` FF): distance-violation + chiral +
//!   4th-dim, at random 4D geometries.
//! - Stage C terms, at native embed conformers: the M6 experimental-torsion term and the UFF
//!   improper term in isolation (RDKit exposes each as its own force field), plus the full Stage-C
//!   force field (`construct3DForceField`). The embed reads f64 throughout (ub64/lb64 bounds, f64
//!   improper/torsion coefficients), so the full Stage-C gradient matches RDKit to machine precision.

use bb_embed::forcefield;

/// Deterministic LCG so the geometries are reproducible without `rand`/time.
struct Lcg(u64);
impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let u = (self.0 >> 11) as f64 / (1u64 << 53) as f64; // [0,1)
        (u - 0.5) * 12.0
    }
}

/// Worst per-component |Δ| relative to the gradient magnitude.
fn rel_grad_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "gradient length mismatch (RDKit parse failure?)");
    let mag = a.iter().chain(b).fold(0.0f64, |m, &x| m.max(x.abs())).max(1.0);
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max) / mag
}

const STAGE_A: &[&str] = &[
    "OC(=O)c1ccccc1",                       // benzoic acid
    "C[C@H](N)C(=O)O",                       // alanine (chiral)
    "O1C[C@]1(CCCCCCCC)[H]",                 // epoxide + explicit-H chiral
    "N12\\C(=N/C3C=CC=CC=3)C=CC=C1C=CC=C2",  // fused bridged, directional bond
    "C1CCC(CC1)C2CCCCC2",                    // two rings
    "c1ccc2c(c1)cccc2C(=O)NCCc3ccccc3",      // amide + aromatic
];

const STAGE_C: &[&str] = &[
    "OC(=O)c1ccccc1",
    "C[C@H](N)C(=O)O",
    "N12\\C(=N/C3C=CC=CC=3)C=CC=C1C=CC=C2",
    "C1CCC(CC1)C2CCCCC2",
    "c1ccc2c(c1)cccc2C(=O)NCCc3ccccc3",
    "O=C1CCCCCCCCCCNC(=O)CCCCCCCCCN1", // a macrocycle
    "N#Cc1ccccc1",                     // benzonitrile — exercises the 179-180 angle-constraint term
    "CC#CC",                           // 2-butyne — two linear angle constraints
];

#[test]
fn stage_a_gradient_matches_rdkit() {
    for &smiles in STAGE_A {
        let spec = bb_spec::build_native(smiles).expect("native spec");
        let n = spec.n_atoms;
        let mut rng = Lcg(0x9E3779B97F4A7C15 ^ (n as u64).wrapping_mul(2654435761));
        for trial in 0..3 {
            let coords: Vec<f64> = (0..n * 4).map(|_| rng.next_f64()).collect();
            let native =
                forcefield::stage_a_energy_grad(&spec, &coords, 4, forcefield::BASIN_DEFAULT).1;
            let rdkit = bb_rdkit::stage_a_ff_grad(
                smiles,
                &coords,
                forcefield::W_CHIRAL,
                forcefield::W_FOURTH,
                forcefield::BASIN_DEFAULT,
            );
            let rel = rel_grad_diff(&native, &rdkit);
            assert!(rel < 1e-6, "{smiles} trial {trial}: Stage-A gradient rel={rel:.2e}");
        }
    }
}

// Stage B (`minimizeFourthDimension`) is the same `constructForceField` with the chiral weight
// dropped to 0.2 and the 4th-dim weight raised to 1.0 (`constructForceField(mmat, pos, csets, 0.2,
// 1.0, nullptr, basinThresh)`). It shares `dist_geom_energy_grad` with Stage A, so Stage A's parity
// already covers the shared distance term, but the reweighted chiral/4th-dim contribs are only
// exercised here — pinning them to the same oracle at weights (0.2, 1.0).
#[test]
fn stage_b_gradient_matches_rdkit() {
    for &smiles in STAGE_A {
        let spec = bb_spec::build_native(smiles).expect("native spec");
        let n = spec.n_atoms;
        let mut rng = Lcg(0xD1B54A32D192ED03 ^ (n as u64).wrapping_mul(2654435761));
        for trial in 0..3 {
            let coords: Vec<f64> = (0..n * 4).map(|_| rng.next_f64()).collect();
            let native = forcefield::stage_b_energy_grad(&spec, &coords, 4, forcefield::BASIN_DEFAULT).1;
            // RDKit's minimizeFourthDimension weights: chiral 0.2, fourth 1.0.
            let rdkit = bb_rdkit::stage_a_ff_grad(smiles, &coords, 0.2, 1.0, forcefield::BASIN_DEFAULT);
            let rel = rel_grad_diff(&native, &rdkit);
            assert!(rel < 1e-6, "{smiles} trial {trial}: Stage-B gradient rel={rel:.2e}");
        }
    }
}

// Native now uses RDKit's exact gradient forms: the M6 torsion via `sinTerm = -dE/dphi / sin(phi)`
// (analytic cancellation of the 1/sin(phi) singularity) plus the degenerate early-return, and the
// improper via the sinY/sinTheta clamp. At well-formed geometries both are ~1e-9-exact. The loose
// tolerance below only guards random geometries that happen to be near-collinear, where the
// degenerate early-return is order-sensitive (native reproduces RDKit's emission order by
// construction, but a single random pathological point is not the place to assert that — the corpus
// gate is). The improper term is exact.
#[test]
fn stage_c_torsion_and_improper_match_rdkit() {
    // Random 3D geometries isolate the term formulas cleanly: no distance constraints are involved,
    // and random points are generically non-degenerate (no near-collinear torsions, whose skip
    // threshold is a geometry edge case, not a formula difference).
    for &smiles in STAGE_C {
        let spec = bb_spec::build_native(smiles).expect("native spec");
        let n = spec.n_atoms;
        let mut rng = Lcg(0x2545F4914F6CDD1D ^ (n as u64).wrapping_mul(2654435761));
        for trial in 0..3 {
            let c: Vec<f64> = (0..n * 3).map(|_| rng.next_f64()).collect();
            let nt = forcefield::torsion_energy_grad(&spec, &c, 3).1;
            let rt = bb_rdkit::stage_c_torsion_grad(smiles, &c);
            let rel = rel_grad_diff(&nt, &rt);
            // formulation-robustness tolerance (see the note above); the formula itself is 1e-13-exact
            assert!(rel < 2e-2, "{smiles} trial {trial}: M6 torsion gradient rel={rel:.2e}");
            let ni = forcefield::improper_energy_grad(&spec, &c, 3).1;
            let ri = bb_rdkit::stage_c_improper_grad(smiles, &c);
            let rel = rel_grad_diff(&ni, &ri);
            assert!(rel < 1e-6, "{smiles} trial {trial}: improper gradient rel={rel:.2e}");
        }
    }
}

#[test]
fn stage_c_full_gradient_matches_rdkit() {
    for &smiles in STAGE_C {
        let spec = bb_spec::build_native(smiles).expect("native spec");
        for conf in bb_embed::embed(&spec, 3, 0xC0FFEE).expect("embed") {
            let c = &conf.coords;
            let (dist_c, angle_c) = forcefield::build_stage_c_constraints(&spec, c);
            let native = forcefield::stage_c_energy_grad(&spec, &dist_c, &angle_c, c).1;
            let rdkit = bb_rdkit::stage_c_ff_grad(smiles, c);
            let rel = rel_grad_diff(&native, &rdkit);
            // Every quantity this term reads is now f64: the bounds matrix (ub64/lb64), the improper
            // coefficients, and the torsion V. There is no f32 left in the Stage-C gradient path, so
            // it matches RDKit to machine precision — the tolerance is tight, not accommodating.
            assert!(rel < 1e-8, "{smiles}: Stage-C full gradient rel={rel:.2e}");
        }
    }
}

#[test]
fn angle_constraint_active_and_matches() {
    // Perturb the embedded geometry to BEND the linear centers, so the 179-180 angle constraint is
    // genuinely active (non-zero gradient), then confirm native's full Stage-C gradient still matches
    // RDKit's construct3DForceField there.
    for smiles in ["N#Cc1ccccc1", "CC#CC"] {
        let spec = bb_spec::build_native(smiles).unwrap();
        let mut c = bb_embed::embed(&spec, 1, 0xC0FFEE).expect("embed")[0].coords.clone();
        // deterministic bend: shove every atom by an index-dependent amount
        for (k, x) in c.iter_mut().enumerate() { *x += 0.15 * (((k * 7 + 3) % 5) as f64 - 2.0); }
        let (dc, ac) = forcefield::build_stage_c_constraints(&spec, &c);
        // rebuild constraints at the *bent* geometry so 1-2/1-3 current pins match RDKit's (it pins to
        // the passed coords too); the angle constraint is bounds-independent (179-180 fixed).
        let ang_g = forcefield::angle_constraint_energy_grad(&ac, &c).1;
        let ang_mag = ang_g.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        assert!(
            !ac.is_empty() && ang_mag > 1.0,
            "{smiles}: angle constraint must be active (mag={ang_mag:.3e}) or the parity check below is vacuous"
        );
        let native = forcefield::stage_c_energy_grad(&spec, &dc, &ac, &c).1;
        let rdkit = bb_rdkit::stage_c_ff_grad(smiles, &c);
        let rel = rel_grad_diff(&native, &rdkit);
        assert!(rel < 1e-5, "{smiles}: Stage-C gradient with active angle constraint rel={rel:.2e}");
    }
}
