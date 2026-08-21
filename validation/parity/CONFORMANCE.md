# Conformance contract

What "faithful reimplementation" means for BetterBuilder, **specified per stage before validating it** —
the thing the field does first and we improvised. RIGOR.md says *how well* a gate must hold; this says
*what relation to the oracle each stage must satisfy, and why*. Every gate cites a row here.

Status of each call is one of: **DECIDED** (ratified, gated), **RECOMMENDED** (proposed, awaiting
ratification), **OPEN** (undecided / unbuilt).

---

## The thesis, and its one honest ceiling

Goal: reproduce, more efficiently, the results Divya's pipeline produces (RDKit ETKDGv3 + her amide
patch → conformer ensembles → charges → db2). The deterministic stages admit exact conformance. The
**stochastic embed core does not** — you cannot prove two RNG-driven pipelines produce "the same"
output. There are exactly two rigorous targets, and they are mutually exclusive:

- **(A) Exact reproduction** — bit-replicate her stochastic core (RDKit Mersenne-Twister draw order +
  her power-iteration eigensolver + her BFGS trajectory + her filters). Yields provable byte-equal
  conformers, but **forces adopting her approximations** — notably the power-iteration eigensolver we
  replaced with an exact one because it is more correct ([[embed-stochastic-layer-decisions]]).
- **(B) Distributional equivalence** — keep our (better) numerics; make a *statistical* claim: native's
  ensemble is indistinguishable from hers within *her own* between-independent-seed variance, with
  confidence bounds. Rigorous statistics, never certainty.

**This choice is the pivot** (see §4). It changes only the embed-core rows; everything else is identical.

---

## 1. Equivalence-relation vocabulary

| Relation | Meaning | Acceptance basis (always *derived*, never picked) |
|---|---|---|
| **byte** | exact format (integers, canonical strings, file bytes) | bit-for-bit |
| **machine-precision** | same formula in f64 | IEEE round-off bound from a *named* cause (f64 same-formula → ~1e-12; a stated cancellation → its condition number) |
| **representation-equivalent** | native uses a different but provably-equivalent encoding | prove the bijection, gate up to it ([[stereo-representation-equivalence]]) |
| **distributional** | stochastic output | indistinguishable from the reference within the reference's *own* between-independent-seed spread; two-sided (coverage + diversity), per-molecule + aggregate, bootstrap CIs |

The oracle is always **the reference executing** — her patched RDKit (build-certified, §3) or her
pipeline output — never our code, never a hand-authored value ([[oracle-is-divyas-pipeline]]).

---

## 2. Per-stage contract

### Setup — `bb-spec` (native) vs her patched RDKit

| Stage | Input domain | Relation | Acceptance | Status |
|---|---|---|---|---|
| atom set / AddHs, formal charge | SMILES | byte | exact atom count & per-atom fields | DECIDED |
| bonds, valence, aromaticity, kekulization, SSSR | SMILES | byte | exact | DECIDED |
| hybridization | SMILES | byte | exact | DECIDED |
| stereo (chiral tags, bond stereo) | SMILES | representation-equivalent | proven bijection to RDKit encoding | DECIDED |
| bounds matrix (12/13/14/15, raw_bounds) | SMILES | machine-precision | f64 same-formula, ~1e-9 (f32-storage floor where it applies) | DECIDED |
| experimental torsions **incl. amide patch** | SMILES | machine-precision | f64; force constants incl. her 80.0/130.9 | DECIDED |
| chiral sets, impropers | SMILES | machine-precision | f64 | DECIDED |

No A/B: setup is deterministic; target is exact / machine-precision. Oracle certified against her
actual build for the FF-relevant path (§3).

### Deterministic per-geometry functions — `bb-embed` vs her patched RDKit

| Stage | Input domain | Relation | Acceptance (derived) | Status |
|---|---|---|---|---|
| Stage-A/B force-field gradient | (mol, 4D geom) | machine-precision | 1e-12 (all-f64 same formula) | DECIDED |
| Stage-C gradient (torsion+improper+constraints) | (mol, 3D geom) | machine-precision | 1e-6 — the V=100 near-planar torsion `1/sinφ` catastrophic-cancellation floor (condition-number-derived), else 1e-12 | DECIDED |
| reject energy (`stage_a_reject_energy`) | (mol, 4D geom) | machine-precision | 1e-12 (RDKit `calcEnergy` identity: native + chiral+4th) | DECIDED |
| 6 post-embed checks | (mol, 3D conformer) | byte | pass/reject bit agreement | DECIDED |
| coord_map bounds tightening | (mol, pins) | machine-precision | 1e-9 (f64 raw_bounds + coordMap) | DECIDED |

No A/B: deterministic given the geometry. **Coverage requirement (Phase 2):** these must be exercised
by divergence-hunting generation over the feature space (ring sizes, amide SMARTS, stereo, near-degenerate
geometries), not only the fixed corpus — the `coord_map` and degenerate-torsion bugs surfaced only on
specific topologies/geometries.

### Stochastic embed core — `bb-embed` init + minimize + accept/reject/retry — **THE A/B DECISION**

Decompose by what each sub-part depends on:

| Sub-part | Depends on | Under (A) exact | Under (B) distributional | Status |
|---|---|---|---|---|
| metric-matrix T construction | distance matrix | same T (identical) | same T (identical) | DECIDED (T is exact either way) |
| eigenstep (T → coords) | the T | port RDKit **power iteration** (magnitude-select, seeded, non-converge) | **keep exact** solver; gated by embeddable-geometry recovery, not RDKit-parity | **DECIDED (B)** |
| distance-matrix draw + neg-eig fill | RNG | match RDKit Mersenne-Twister draw order | **any good RNG**; validated only at the ensemble | **DECIDED (B)** |
| minimizer | init | match BFGS + backtracking trajectory | **our L-BFGS** (proven energy-equivalent from a good init) | **DECIDED (B)** |
| reject filters: per-atom energy (0.05), eigenvalue `numZeroFail`, retry cap | geometry / eigenvalues | port faithfully | **port faithfully** (they shape the accepted *distribution*) | DECIDED — needed under both |
| conformer output | all the above | byte-equal to her seeded conformer | **distributional match** to her ensemble, CALIBRATED to her own run-to-run: native's per-molecule pass rate on the strict recovery+diversity gate ≥ her OWN independent run's (base2). Fairness control: native 29/96 (30%) vs her base2 17/94 (18%) on the IDENTICAL gate → native is at least as self-consistent with her ensemble as she is with herself = 'within her run-to-run variability.' The absolute ≥95% is unachievable even by her own pipeline (18%) at this strictness, so the relation is the calibrated one. | **DECIDED (B) — ACCEPT, calibrated gate GREEN (native ≥ her run-to-run)** |

The reject filters are **not** an A/B axis — they are faithful behavior we currently lack and must
port either way, because they determine which draws survive into the ensemble.

### Output — `bb-output` / `bb-db2` vs her pipeline output

| Stage | Input | Relation | Acceptance | Status |
|---|---|---|---|---|
| db2 | fixed conformer | byte | exact vs her real files (3287) | DECIDED ([[db2-output-gap-map]]) |
| mol2 / SYBYL types / strain | fixed conformer | byte | exact vs her output | DECIDED / partial |
| AMSOL charges | conformer from her **unseeded** solvation embed | **distributional / tolerance** | RESOLVED by investigation: her solvation embed IS unseeded (`build_ligands.py` L112 `EmbedMolecule(mol, ETKDGv3())`, no randomSeed → -1 → non-deterministic) AND AMSOL charges are geometry-SENSITIVE (measured per-atom CM2 spread across conformers: max 0.21e, mean 0.047e). So her OWN charges vary run-to-run by ~0.2e — byte-parity is impossible even for her pipeline against itself; the only coherent relation is tolerance-based within her ~0.2e envelope. Same distributional situation as the embed core. | **DECIDED — tolerance (~0.2e), byte-parity impossible by her own non-determinism** |

Output parity is *conditional on the input conformer*; given a fixed conformer it is deterministic
(byte). The conformer itself comes from the stochastic tier, so end-to-end output parity inherits the
A/B outcome.

---

## 3. Oracle provenance (must be complete for the contract to hold)

Every deterministic golden is only as faithful as "our oracle == her actual build."

- **Established (standing `gates.sh certify` tier, `validation/certify/`):** her production container's
  RDKit is `2026.09.1pre` (= ours); her build's compiled `libRDKitForceFieldHelpers` carries her amide
  patch (`130.9`); her ETKDGv3 config is stock but for `randomSeed` + `useExpTorsionAnglePrefs`
  (already default). Oracle dumpers compiled against **her actual `/soft/rdkit_libs`** and diffed vs
  our host-bridge goldens are **byte-identical (`0.00e+00`) over the full 100-molecule corpus** for
  **all six deterministic oracles**: Stage-A (base FF), Stage-B (`constructForceField` at 0.2/1.0,
  same dumper), Stage-C (`getExperimentalTorsions` incl. the amide patch — the amide-patch proof),
  the six post-embed **checks**, the smoothed **bounds** matrix, and the **coord_map**
  (`adjustBoundsMatFromCoordMap`) tightening. Our from-source build == her production build, proven
  across every deterministic component.
- **Phase 1 COMPLETE.** Every deterministic golden is now Divya-faithful by proof, not by patch
  reconstruction. The standing `gates.sh certify` tier re-proves it on demand.

---

## 4. The pivotal decision — RATIFIED: (B) distributional

**Decision: (B) distributional equivalence, with the deterministic reject filters ported.** Ratified;
the embed-core and output-conformer rows are now DECIDED against target (B). We keep the exact
eigensolver and our minimizer, validate the embed only at the *ensemble* level (never
conformer-for-conformer), and port the reject/retry filters so the accepted distribution matches hers.

Rationale:
1. **(A) forces reproducing a defect.** It requires downgrading the exact eigensolver to RDKit's
   seeded power iteration, which we established is genuinely less correct (magnitude-selects
   non-Euclidean modes, fails to converge). That violates strict-improvement / never-reproduce-her-defects
   ([[embed-stochastic-layer-decisions]], [[no-shortcuts-ever]]).
2. **(B) is achievable on the evidence.** The eigensolver difference lands *within* intrinsic sampling
   spread geometrically, her own 0.05/atom energy reject culls the divergent bad-init outcomes, and our
   minimizer is energy-equivalent to hers from a good init — so a distributional match is plausible
   once the filters are ported.
3. **"Her results" most usefully means her distribution** — the ensemble's coverage/diversity that
   downstream consumes — not the exact coordinates of conformer #37.

**What ratifying (B) commits us to:** port the reject/retry filters; build the end-to-end distributional
gate (native ensemble vs her *actual* seeded ensembles, acceptance = her own independent-seed spread,
symmetry-corrected RMSD + TFD + strain, per-molecule + bootstrap CIs); keep the exact solver and our
minimizer, validated *only* at the ensemble, not conformer-for-conformer.

**What ratifying (A) instead would commit us to:** replicate RDKit's Mersenne-Twister draw order, its
power-iteration eigensolver, and its BFGS trajectory bit-for-bit, accepting the numerics downgrade, in
exchange for byte-equal conformers.

**Overrule this if** the true downstream requirement is her *exact* conformers (then (A)).

---

## 5. What changes once ratified

- Phase 1 — complete oracle certification against her container (§3) → every deterministic golden proven.
- Phase 2 — replace corpus-sampling with divergence-hunting generation for every deterministic row.
- Phase 3 — build the embed-core gate for the ratified target (A: bit-match; B: distributional).
- Phase 4 — close the output tier to byte-parity, resolving the unseeded-solvation-embed question.

Until §4 is ratified, the embed-core and output-conformer rows are **OPEN by definition** — not for
lack of effort, but because the target is unspecified.

---

## 6. Measured status (seeds_100 corpus)

- **Phase 1 — COMPLETE.** All six deterministic oracles byte-identical (`0.00e+00`) to her *actual
  container* build over 100 molecules (`validation/certify/`, `gates.sh certify`). Answers "are we
  considering her patches properly": our from-source build == her production build.
- **Phase 3 — embed-core distributional gate, RUNNING GREEN on RMSD.** `check.py stage3` over 96 clean
  ensembles (her verbatim recipe at seed base 210185 vs native `embed_recipe`), each metric normalised
  by her OWN independent-half spread:
  - RMSD NN **1.07×** [1.04, 1.09], within 1.08×, between 1.06× — native reproduces her ensemble to
    within ~7% of an independent run of *her own* pipeline. Recovery curve overlays hers.
  - TFD NN **1.44×** [1.33, 1.58] — elevated (ring-torsion metric); its acceptance is decided by the
    derived band, not asserted.
  - Aggregate-CI gate PASS; per-molecule 75/96 under the *placeholder* ±15% band (most fails marginal).
- **Derived band (STEP 4) — mechanism built, second-base run in progress.** `check.py band` reruns her
  verbatim recipe at an independent seed base (918273) and measures the SAME ratios base-1-vs-base-2 →
  her run-to-run variability → `derived_band.json`, which stage3 auto-consumes in place of the ±15%
  placeholder. This is the "measured band, not placeholder" the target requires.
- **Dims 3+4 (improvement / speed) — MEASURED.** Native median 0.336 s/molecule (RDKit-free) → ~119×
  over her ~40 s; staller census 4/100 (seed-specific infeasible coordMap), native emits a full
  200-conf ensemble on all 4 where she emits 0. Confirmed seed-specificity: a base-1 staller
  (CSLB00000xi1qR) embeds cleanly (200 confs, 14.8 s) at base 918273.
- **Phase 2 — divergence-hunting BUILT + GREEN (deterministic layers).** `validation/hunt/`
  (`gates.sh hunt`): a structure-aware generator mutates her real `seeds_5000` workload into adversarial
  macrocycles (heteroatom swap / substituent add / terminal delete / bond-order flip, valence-valid,
  ring ≥ 9, fixed-seed → deterministic), with **measured coverage** (ring 13–22, heteroatoms 5–21,
  atoms 28–51). **Spec layer: 6832 generated mutants diffed native vs the RDKit bridge → 0 divergences
  of every kind** (0 parse, 0 hard = elements+AddHs-order/bonds/recipe/torsion+improper+angle topology
  incl. amide patch, 0 numeric, 0 stereo — the 7 `stereo_double_bonds` flags proven pure reference-atom
  convention by a representation-invariant canonical form, not asserted). **Embed layer: 500 mutants
  through the vetted `corpus_parity` gate → all machine-precision** (Stage-A 3.7e-16, B 4.2e-16,
  C 1.0e-10, checks 0/500, coord_map 0.0, reject 4.0e-15). So "no diffs" now means "searched hard,"
  not "the sample passed" — closing criterion 5's generative-coverage gap for the deterministic layers.
- **Phase 4 — output tier OPEN.** `bb-db2` writer byte-exact on the *format* (9 real files round-trip),
  but the full path native-conformers → mol2 → AMSOL charges → db2 (the mol2db2 hierarchy) is unbuilt;
  AMSOL charge determinism (unseeded solvation embed) still to resolve. This is the largest remaining gap.
