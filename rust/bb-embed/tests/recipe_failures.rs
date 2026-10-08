//! The core-pin recipe when embeds fail (GitHub issue #1: a bundle crashed with an index-out-of-
//! bounds panic in `embed_recipe` when a core seed embedded no conformer).
//!
//! A caller who asks for N conformers gets N accepted conformers. Every embed gets RDKit's attempt
//! budget, and only conformers that pass the acceptance checks are returned. The recipe returns its
//! blocks, each built on one core; a block that cannot be completed is rebuilt around the next unused
//! core, and a failed independent conformer is drawn again. When every allowed replacement or redraw
//! also fails, or the spec lacks an input the recipe needs, the result is an error, never a panic or a
//! short ensemble.

use bb_core::{MoleculeSpec, StereoDoubleBond};
use bb_embed::{Block, Conformer, RecipeError, Shortfall};

/// The molecule from issue #1 (CSLB0000477NCT).
const ISSUE_1_SMILES: &str =
    "COC[C@H]1C[C@H]2CNC(=O)CC(C)C=CCCC(NC(=O)NC3=NN(CCC(=O)OC(C)(C)C)C=N3)C(=O)N2C1";
/// mc0405 from the 1000-macrocycle test set (`macrocycles_21-50HA_1000.smi`).
const MC0405_SMILES: &str = r"O=C1/C=C(\CF)COC(CF)CCCOC(=O)C[C@H](CS)O1";
/// mc0939 from the same test set.
const MC0939_SMILES: &str = "O=C1CC[C@@H](CS)OC(=O)CC[C@H](O)COC(CCl)CO1";
/// `bb-build`'s default seed, under which issue #1's molecule panicked.
const BB_BUILD_SEED: u64 = 0x0BAD_C0DE;

/// Three atoms whose bounds fix d01 = d12 = 1 Å and d02 = 10 Å. No geometry in any dimension
/// meets them: if d01, d12 ≤ 2 then d02 ≤ 4 by the triangle inequality, and the d02 term of RDKit's
/// distance-violation energy alone is ≥ (2·10²/(10²+4²) − 1)² ≈ 0.524; otherwise a d01 or d12 term
/// is ≥ (2² − 1)² = 9. Both exceed RDKit's `firstMinimization` reject threshold of
/// 0.05 × 3 atoms = 0.15, so every embed attempt is rejected and no conformer can exist.
fn infeasible_spec(core_seeds: u32, sidechain_confs: u32) -> MoleculeSpec {
    let n = 3;
    let target = [[0.0, 1.0, 10.0], [1.0, 0.0, 1.0], [10.0, 1.0, 0.0]];
    let bounds_f64: Vec<f64> = (0..n * n).map(|k| target[k / n][k % n]).collect();
    MoleculeSpec {
        n_atoms: n,
        dim: 4,
        bounds: bounds_f64.iter().map(|&x| x as f32).collect(),
        raw_bounds: bounds_f64.iter().map(|&x| x as f32).collect(),
        raw_bounds_f64: bounds_f64.clone(),
        bounds_f64,
        pin_atoms: vec![0, 1, 2],
        core_seeds,
        sidechain_confs,
        ..Default::default()
    }
}

/// Three atoms, 0 and 1 pinned, whose sidechain embeds succeed or fail depending on the core. The
/// core reads the smoothed bounds, every distance in [1, 4] Å, so each core embeds with its own d01.
/// The sidechains start from the raw bounds, as RDKit's coordMap path does: d02 = d12 = 1 Å.
/// - Core d01 ≤ 1.8 Å: atom 2 at 1 Å from two atoms ≤ 2 Å apart is realisable, so the sidechains
///   embed.
/// - Core d01 ≥ 2.4 Å: d01 > d02 + d12, so RDKit cannot triangle-smooth and embeds nothing
///   (`setupInitialBoundsMatrix`). A pin-free fallback to the raw bounds (d01 ≤ 1.2 Å) leaves the
///   pinned pair's violation alone at ≥ (2.4²/1.2² − 1)² = 9 > 0.15, RDKit's reject threshold for
///   3 atoms. Either way that core yields no sidechain conformer.
fn gap_spec(core_seeds: u32, sidechain_confs: u32) -> MoleculeSpec {
    let n = 3;
    // Upper triangle = upper bound, lower triangle = lower bound.
    let smoothed = [[0.0, 4.0, 4.0], [1.0, 0.0, 4.0], [1.0, 1.0, 0.0]];
    let raw = [[0.0, 1.2, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]];
    let flat = |m: [[f64; 3]; 3]| -> Vec<f64> { (0..n * n).map(|k| m[k / n][k % n]).collect() };
    let (bounds_f64, raw_bounds_f64) = (flat(smoothed), flat(raw));
    MoleculeSpec {
        n_atoms: n,
        dim: 4,
        bounds: bounds_f64.iter().map(|&x| x as f32).collect(),
        raw_bounds: raw_bounds_f64.iter().map(|&x| x as f32).collect(),
        bounds_f64,
        raw_bounds_f64,
        pin_atoms: vec![0, 1],
        core_seeds,
        sidechain_confs,
        ..Default::default()
    }
}

/// Six disconnected units a–b=c–d, each with its double bond required trans by a stereo double bond,
/// while its 1-4 bound, [2.8, 4.0] Å, admits both cis (≈3.0 Å) and trans (≈3.7 Å) for these bond lengths
/// and 1-3 distances, and no torsion term prefers either. An attempt is accepted only if every unit lands
/// trans, so a draw of 10 × 24 attempts sometimes fails and sometimes succeeds: embeds that need a redraw.
fn trans_units_spec() -> MoleculeSpec {
    let (k, n) = (6, 24);
    let mut b = vec![0.0f64; n * n];
    let mut set = |i: usize, j: usize, l: f64, u: f64| {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        b[lo * n + hi] = u;
        b[hi * n + lo] = l;
    };
    for i in 0..n {
        for j in (i + 1)..n {
            set(i, j, 3.0, 50.0);
        }
    }
    let mut stereo_double_bonds = Vec::new();
    for u in 0..k {
        let (a, bb, c, d) = (4 * u, 4 * u + 1, 4 * u + 2, 4 * u + 3);
        set(a, bb, 1.5, 1.5);
        set(bb, c, 1.34, 1.34);
        set(c, d, 1.5, 1.5);
        set(a, c, 2.48, 2.48);
        set(bb, d, 2.48, 2.48);
        set(a, d, 2.8, 4.0);
        stereo_double_bonds.push(StereoDoubleBond {
            atoms: [a as u32, bb as u32, c as u32, d as u32],
            sign: 1,
        });
    }
    MoleculeSpec {
        n_atoms: n,
        dim: 4,
        bounds: b.iter().map(|&x| x as f32).collect(),
        bounds_f64: b,
        stereo_double_bonds,
        ..Default::default()
    }
}

fn core(spec: &MoleculeSpec, base: u64, j: usize) -> Conformer {
    bb_embed::recipe_core(spec, base, j).unwrap_or_else(|| panic!("core {j} embeds"))
}

fn d01(c: &Conformer) -> f64 {
    let (a, b) = (c.atom(0), c.atom(1));
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// The recipe's block contract: block `s` is built on `cores[s]` and holds `spec.sidechain_confs`
/// complete conformers, each with `spec.pin_atoms` exactly at that core's coordinates (they are held
/// fixed, not recomputed) and each passing RDKit's acceptance checks against the bounds it was embedded
/// to: the raw bounds tightened to the core's pins (RDKit's coordMap path).
fn assert_blocks(spec: &MoleculeSpec, base: u64, blocks: &[Block], cores: &[usize]) {
    let n = spec.n_atoms;
    let built: Vec<usize> = blocks.iter().map(|b| b.core).collect();
    assert_eq!(built, cores, "the blocks' cores");
    for b in blocks {
        assert_eq!(
            b.conformers.len(),
            spec.sidechain_confs as usize,
            "block on core {}",
            b.core
        );
        let core = core(spec, base, b.core);
        let mut pinned = vec![None; n];
        for &a in &spec.pin_atoms {
            pinned[a as usize] = Some(core.atom(a as usize));
        }
        let sc_bounds = bb_embed::bounds::coord_map_bounds_f64(&spec.raw_bounds_f64, n, &pinned);
        let sc_spec = spec.with_seed_bounds(sc_bounds);
        for (i, c) in b.conformers.iter().enumerate() {
            assert_eq!(
                c.coords.len(),
                n * 3,
                "core {} conformer {i} is incomplete",
                b.core
            );
            assert!(
                c.coords.iter().all(|x| x.is_finite()),
                "core {} conformer {i}",
                b.core
            );
            for &a in &spec.pin_atoms {
                assert_eq!(
                    c.atom(a as usize),
                    core.atom(a as usize),
                    "core {} conformer {i}, pin {a}",
                    b.core
                );
            }
            assert!(
                bb_embed::checks::passes_checks(&sc_spec, &c.coords, &pinned),
                "core {} conformer {i} fails the acceptance checks",
                b.core
            );
        }
    }
}

#[test]
fn recipe_rejects_a_spec_without_recipe_counts() {
    for (core_seeds, sidechain_confs) in [(0, 3), (2, 0), (0, 0)] {
        match bb_embed::embed_recipe(&gap_spec(core_seeds, sidechain_confs), BB_BUILD_SEED) {
            Err(RecipeError::NoRecipeCounts {
                core_seeds: c,
                sidechain_confs: s,
            }) => {
                assert_eq!((c, s), (core_seeds, sidechain_confs));
            }
            other => {
                panic!("({core_seeds}, {sidechain_confs}): expected NoRecipeCounts, got {other:?}")
            }
        }
    }
}

#[test]
fn recipe_rejects_a_spec_without_raw_bounds() {
    let mut spec = gap_spec(2, 2);
    spec.raw_bounds.clear();
    spec.raw_bounds_f64.clear();
    assert!(matches!(
        bb_embed::embed_recipe(&spec, BB_BUILD_SEED),
        Err(RecipeError::NoRawBounds)
    ));
}

#[test]
fn recipe_reports_a_shortfall_when_no_core_can_embed() {
    match bb_embed::embed_recipe(&infeasible_spec(2, 3), BB_BUILD_SEED) {
        Err(RecipeError::Shortfall(Shortfall {
            requested,
            embedded,
        })) => {
            assert_eq!((requested, embedded), (6, 0));
        }
        other => panic!("expected a shortfall, got {other:?}"),
    }
}

#[test]
fn embed_returns_every_requested_conformer() {
    let confs =
        bb_embed::embed(&gap_spec(0, 0), 7, BB_BUILD_SEED).expect("a realisable spec embeds");
    assert_eq!(confs.len(), 7);
    assert!(confs
        .iter()
        .all(|c| c.coords.len() == 9 && c.coords.iter().all(|x| x.is_finite())));
}

#[test]
fn embed_reports_a_shortfall_when_nothing_can_embed() {
    let shortfall = bb_embed::embed(&infeasible_spec(0, 0), 3, BB_BUILD_SEED)
        .expect_err("no conformer of this spec exists");
    assert_eq!((shortfall.requested, shortfall.embedded), (3, 0));
}

/// A conformer whose draw fails is drawn again. `recipe_core(spec, s, 0)` is `embed`'s first draw for
/// seed `s` (the same seed and init), so with that draw failing, a conformer can only come from a redraw.
#[test]
fn embed_redraws_a_conformer_whose_draw_fails() {
    let (spec, seed) = (trans_units_spec(), 2);
    assert!(
        bb_embed::recipe_core(&spec, seed, 0).is_none(),
        "the first draw must fail"
    );
    let confs = bb_embed::embed(&spec, 1, seed).expect("the redraw embeds");
    assert_eq!(confs.len(), 1);
    assert!(bb_embed::checks::passes_checks(
        &spec,
        &confs[0].coords,
        &[]
    ));
}

/// Redrawing stops at the safeguard: 2 draws for one conformer. For this seed draws 1 and 2 fail and draw
/// 3 would succeed (measured with this engine, so the test pins the safeguard rather than proving it), so
/// a third draw would complete the request and the safeguard must report a shortfall instead.
#[test]
fn embed_stops_redrawing_at_the_safeguard() {
    assert_eq!(
        bb_embed::MAX_TRIES_PER_REQUEST,
        2,
        "the seed below is chosen for 2 draws per conformer"
    );
    let (spec, seed) = (trans_units_spec(), 1);
    assert!(
        bb_embed::recipe_core(&spec, seed, 0).is_none(),
        "the first draw must fail"
    );
    let shortfall = bb_embed::embed(&spec, 1, seed).expect_err("both allowed draws fail");
    assert_eq!((shortfall.requested, shortfall.embedded), (1, 0));
}

/// A core that cannot hold its sidechains is replaced by the next unused core, the k-th such block
/// taking core `core_seeds + k`, and every block comes out complete and in block order.
#[test]
fn recipe_replaces_a_core_that_cannot_hold_its_sidechains() {
    let (spec, base) = (gap_spec(6, 2), BB_BUILD_SEED + 106);
    // Preconditions: each core used falls clearly on one side of the gap; cores 2 and 5 cannot hold
    // their sidechains, and their replacements, cores 6 and 7, can.
    let held: Vec<bool> = (0..8)
        .map(|j| {
            let d = d01(&core(&spec, base, j));
            assert!(
                d <= 1.8 || d >= 2.4,
                "core {j}: d01 = {d} Å is inside (1.8, 2.4)"
            );
            d <= 1.8
        })
        .collect();
    assert_eq!(held, [true, true, false, true, true, false, true, true]);

    let blocks = bb_embed::embed_recipe(&spec, base).expect("the recipe completes");
    assert_blocks(&spec, base, &blocks, &[0, 1, 6, 3, 4, 7]);
}

/// Replacement stops after `core_seeds` replacement cores. With one core: core 0 and its one allowed
/// replacement, core 1, cannot hold their sidechains, so the recipe reports a shortfall even though
/// core 2, which it is not allowed to try, could.
#[test]
fn recipe_stops_after_core_seeds_replacements() {
    let (spec, base) = (gap_spec(1, 1), BB_BUILD_SEED + 30);
    let d: Vec<f64> = (0..3).map(|j| d01(&core(&spec, base, j))).collect();
    assert!(
        d[0] >= 2.4 && d[1] >= 2.4,
        "cores 0 and 1 must be unable to hold: d01 = {d:?}"
    );
    assert!(d[2] <= 1.8, "core 2 must be able to hold: d01 = {d:?}");
    match bb_embed::embed_recipe(&spec, base) {
        Err(RecipeError::Shortfall(Shortfall {
            requested,
            embedded,
        })) => {
            assert_eq!((requested, embedded), (1, 0));
        }
        other => panic!("expected a shortfall, got {other:?}"),
    }
}

/// trans-Cycloheptene as written, `C1CC/C=C/CC1`: RDKit ignores E/Z on a double bond in a ring smaller
/// than 8, so the oracle embeds it (RDKit `EmbedMolecule`: 50 of 50 seeds, measured 2026-10-06). Before
/// perception applied that rule, every attempt failed the trans check and the recipe reported a
/// shortfall; now the spec carries no stereo double bond and the recipe completes.
#[test]
fn a_small_ring_stereo_marker_is_ignored_as_rdkit_does() {
    let spec = bb_spec::build_native("C1CC/C=C/CC1").expect("spec");
    assert!(
        spec.stereo_double_bonds.is_empty(),
        "no E/Z on a 7-ring double bond"
    );
    let blocks = bb_embed::embed_recipe(&spec, BB_BUILD_SEED).expect("the recipe completes");
    let cores: Vec<usize> = blocks.iter().map(|b| b.core).collect();
    assert_eq!(cores.len(), spec.core_seeds as usize);
    assert_blocks(&spec, BB_BUILD_SEED, &blocks, &cores);
}

/// Issue #1's molecule at `bb-build`'s seed. Its 10th core exhausted the former 24-attempt budget,
/// and the recipe panicked; with RDKit's budget of 10 × atoms every core and sidechain embeds, so
/// every block is built on its own core. That every core embeds is this engine's measured behaviour
/// (the 10th core needs 37 attempts), not an independent result: the test pins it against regression,
/// while the block contract it asserts comes from the recipe.
#[test]
fn issue_1_molecule_embeds_on_its_own_cores() {
    let spec = bb_spec::build_native(ISSUE_1_SMILES).expect("issue #1 SMILES builds a spec");
    let blocks = bb_embed::embed_recipe(&spec, BB_BUILD_SEED).expect("the recipe completes");
    let own: Vec<usize> = (0..spec.core_seeds as usize).collect();
    assert_blocks(&spec, BB_BUILD_SEED, &blocks, &own);
}

/// mc0405 at `bb-build`'s seed: its sidechain conformer (core 0, k = 15) needs 228 attempts, more
/// than a 24-attempt cap or 10 + 42 atoms = 52 attempts, and fewer than RDKit's 10 × 42 = 420. A
/// sidechain's RNG depends only on (base, core, k), so the recipe cut to 1 core × 16 sidechains embeds
/// it unchanged. Under RDKit's budget it succeeds, so core 0 keeps its block; a smaller budget would
/// fail it and replace the core. The 228 attempts are this engine's measured behaviour, not an
/// independent result: the test pins the budget against regression rather than proving it.
#[test]
fn a_slow_sidechain_within_rdkits_budget_keeps_its_core() {
    let mut spec = bb_spec::build_native(MC0405_SMILES).expect("mc0405 SMILES builds a spec");
    (spec.core_seeds, spec.sidechain_confs) = (1, 16);
    let blocks = bb_embed::embed_recipe(&spec, BB_BUILD_SEED).expect("the recipe completes");
    assert_blocks(&spec, BB_BUILD_SEED, &blocks, &[0]);
}

/// mc0939 at `bb-build`'s seed: no sidechain embeds around core 5. Real RDKit (2026.09.1pre with
/// the amide patch, `EmbedMultipleConfs` with core 5's pins as its coordMap, measured 2026-10-05)
/// embedded 0 of 20 there and 20 of 20 around each of cores 0–4 and 6–9. Cut to 6 cores × 4
/// sidechains, the recipe must replace core 5 with core 6, and every conformer it returns must pass
/// the acceptance checks; no rejected attempt is passed off as a conformer.
#[test]
fn a_core_rdkit_cannot_reembed_around_is_replaced() {
    let mut spec = bb_spec::build_native(MC0939_SMILES).expect("mc0939 SMILES builds a spec");
    (spec.core_seeds, spec.sidechain_confs) = (6, 4);
    let blocks = bb_embed::embed_recipe(&spec, BB_BUILD_SEED).expect("the recipe completes");
    assert_blocks(&spec, BB_BUILD_SEED, &blocks, &[0, 1, 2, 3, 4, 6]);
}
