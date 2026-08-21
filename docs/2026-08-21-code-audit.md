# BetterBuilder — code audit

Architectural-debt / performance sweep, 2026-08-21. Covers the whole workspace: the ~10 Rust crates
and the `validation/` layer. Method: three parallel deep reviews (embed hot path · output pipeline ·
validation layer) plus a `cargo clippy` pass. This is a **work plan** — nothing here has been applied
yet, but all of it will be: on this project, organization, single-sourcing, architecture, and optimality
are standards, not optional cleanups. The tiers below sequence the work; they never gate whether it
happens.

Line references are `path:line` as of this date. Each finding is tagged **behaviour-preserving** (the
fix cannot change output bytes — safe under the bit-exact parity posture) or **needs re-validation**
(the fix reassociates floating-point or otherwise perturbs the numerics that drive the minimiser, so it
must pass `corpus_parity` + `stage3` before landing).

Context that shapes the priorities: much of the duplication below was introduced during a single long
build-out of the output pipeline and the criterion-2 gate.

## Correctness defects (verified against the code)

These are bugs. They are listed here, separately and without severity ranking — a bug is not "minor,"
and impact-tier is orthogonal to whether something is incorrect. Every one is in the validation /
benchmark tooling; the shipping product (`bb-embed`'s math, the `bb-output`/`bb-db2` writers) has **no
identified correctness bug** — its findings below are all performance/duplication with the current output
verified correct. "Latent" means the defect is not exercised by the current corpus but is wrong for
inputs the code is written to handle; it is still a bug.

- **C1 — median computed wrong for even-length lists (ACTIVE).** `speedup_report.py:58` (`s[len(s)//2]`)
  and `time_native_recipe.py:56` (`ts[len(ts)//2]`) take the upper of the two middle elements, not their
  average. Every "median" these scripts print for an even-count set is wrong. (Independent of the
  threading issue; the corrected 4.7× used `statistics.median` and is right.) Fix: `statistics.median`.
- **C2 — `mol_id` reimplementations are not equivalent to the canonical (LATENT).** The copies
  (`excess_validity.py:54`, `diagnose_tfd.py:103`, `strain_compare.py:83`, `bounds_strain.py:48`) do only
  `name.split("_")[0]` and omit the canonical `check.mol_id`'s `.<prot_id>` suffix strip. On a name like
  `mol.0` they return `mol.0` where the canonical returns `mol` → wrong molecule key → the molecule is
  silently skipped or mismatched. They coincide only because `seeds_100` names carry no such suffix —
  which the canonical's docstring says "her records" do. Fix: import `check.mol_id`.
- **C3 — criterion-2 `DIV_MIN` fallback disagrees across copies (LATENT).** With `derived_band.json`
  absent, `check.py:376` uses `0.85` while `criterion2_gate.py:31` and `fairness_control.py:30` use
  `0.99` → the same gate can return different verdicts by entry point, and silently uses an unauthoritative
  threshold instead of erroring. Not triggered in any reported run (the band file existed). Fix: one
  loader, one fallback (or error). See finding 1.
- **C4 — fixed-path temp-spec files race under concurrent runs (LATENT).** `bounds_strain.py:43`
  (`_spec_tmp.json`) and `excess_validity.py:70` (`_ev_spec.json`) write a spec to a fixed repo path per
  molecule; two concurrent instances clobber it and one scores conformers against the wrong molecule's
  spec. (`_ev_spec.json` is also never unlinked.) Fix: `tempfile`. See finding 12.

The debt/perf findings below are ranked by engineering impact; that ranking says nothing about
correctness and none of the tiered items are bugs (finding 2, once flagged, was verified NOT to be one).

---

## Summary

Correctness defects are the section above (C1–C4); they are not tiered. The tiers below are debt/perf
only — none of them is a bug — ranked by engineering impact:

| Tier | Theme | Count | Nature |
|---|---|---|---|
| 1 | SSoT hazard | 2 | Finding 1 also underlies correctness defect C3; finding 2 was verified NOT a bug |
| 2 | Hot-path performance | 6 | 5 behaviour-preserving, 1 needs re-validation |
| 3 | DRY / SSoT / architecture | ~9 | Maintainability; findings 11/12 also underlie correctness defects C2/C4 |
| 4 | Small cleanups | ~4 | Genuine cleanups (no bugs — the median defect lives in C1, not here) |

`clippy` on the RDKit-free crates is essentially clean (~10 trivial style warnings: loop-index-instead-
of-iterator, one "slice will do" allocation, a redundant cast). The debt is architectural, not
mechanical — the codebase is well-kept at the lint level.

---

## Tier 1 — correctness / SSoT hazards

### 1. Criterion-2 gate has three copies with a divergent fallback constant
`validation/parity/check.py:376,408,410-413` · `criterion2_gate.py:23,29-31,34-37` · `fairness_control.py:20,26-31,33-37`

`REC_TOL`, the `failed(row, div_min)` body, and the band loader are triplicated. The `failed()` logic is
byte-identical, but the **`DIV_MIN` fallback differs**: with `derived_band.json` absent, `check.py`
falls back to the provisional `0.85`, while both gate scripts fall back to `0.99`. The same "ratified
gate" would therefore apply a different diversity floor depending on which entry point runs it.

**This is correctness defect C3.** It did not affect reported numbers — `derived_band.json` existed in
every run, so all three used the derived `0.99` — but the code is incorrect by construction: a run
without the band file passes under one runner and fails under another. **Fix:** promote `REC_TOL`,
`failed(row, div_min)`, and one `load_band()` (one fallback, or an error) into `check.py`; both drivers
import them.

### 2. Gate and production differ on strain rounding — duplication only, NOT a correctness bug
`rust/bb-output/src/bin/emit_inputs.rs:79-81` vs `rust/bb-output/src/assemble.rs:137-139`

`emit_inputs.rs` (the writer **gate**) rounds per-conformer strain to `%.3f` before the db2 build;
`build_blocks_db2` (**production**) slices strain without that rounding. This was flagged as a possible
correctness blind spot (the gate validating a different artifact than ships). **On verification it is
not a bug.** The S-record formats strain as `{:+11.3}` (`bb-db2/src/write.rs:117`), i.e. `%.3f`, and
pre-rounding to `%.3f` is idempotent under `%.3f` formatting: the gate emits `%.3f(round(x,3)) =
round(x,3)`, production emits `%.3f(x) = round(x,3)` — **byte-identical**. So production and the gate
produce the same db2, and the 100/100 corpus writer byte-identity result *does* cover production. The
gate's extra `%.3f` round exists only to stop the `%.6f` mol2-comment round-trip flipping the parse on a
half-way boundary; it changes nothing on the db2 side.

What remains is real but is **debt, not incorrectness**: the two code paths are a copy-paste (see finding
9) and could drift into a genuine divergence later. **Fix:** single-source the per-molecule pipeline so
they cannot. Behaviour-preserving.

---

## Tier 2 — hot-path performance

The cost centre is the FF eval loop: `embed_recipe` → parallel `embed_one` → `minimize_stage_{a,b,c}` →
argmin L-BFGS calling `*_energy_grad` many times per conformer, ×200 conformers ×N molecules. The inner
math is already well-structured (energy+gradient fused, Stage-C constraints pre-built once per attempt,
eval cache memoised, `grad_tol()` `OnceLock`-memoised so `BB_GRAD_TOL` is *not* re-read per call).
Remaining wins:

### 3. Stage-A/B rebuild the distance-pair list on every FF eval — **the biggest single win**
`rust/bb-embed/src/forcefield.rs:134-172` (`accum_dist_term`) · **behaviour-preserving**

Every eval scans all n²/2 pairs, re-reads `ub64`/`lb64` (each with an `i<j` swap), and re-applies the
`ub-lb > basin` filter — but the surviving pair set and its `lb2`/`ub2` are **constant across the whole
minimisation**. Stage C already pre-builds its contrib list; A/B don't, and A/B run 400+200 iterations.
**Fix:** build `Vec<(i,j,lb2,ub2)>` of basin-surviving pairs once in stage setup, iterate that each eval
(also skips the majority of pairs on sparse molecules instead of branching over them).

### 4. Redundant `sqrt` in that same O(n²) loop
`rust/bb-embed/src/forcefield.rs:149,154,161` · **behaviour-preserving**

`d2.sqrt()` is computed in the `> ub2` and `< lb2` branches and again at line 161. Compute
`let d = d2.sqrt()` once when the pair violates and reuse.

### 5. `build_blocks_db2` deep-clones the whole `TypedMol` per block
`rust/bb-output/src/assemble.rs:121,133` (mirrored `emit_inputs.rs:71`) · **behaviour-preserving**

Each block `mol.clone()`s every conformer + all topology, then immediately overwrites `atom_xyz` with a
`[lo..hi]` slice; there is a second full clone at `:121` to round coords. O(nblocks·nconf·natoms) copies
where O(nconf·natoms) suffices. **Fix:** round coords once into an owned mol, then per block clone only
the topology and move in the `atom_xyz`/strain slices — or have `hierarchy::build` accept slices so no
per-block `TypedMol` is materialised.

### 6. `embed_recipe` clones the entire `MoleculeSpec` per core seed
`rust/bb-embed/src/lib.rs:196` · **behaviour-preserving**

Clones four n² bounds matrices (`bounds`, `bounds_f64`, `raw_bounds`, `raw_bounds_f64`) plus all
torsion/improper/angle/bond arrays, 10–20× per molecule, only to overwrite `bounds`/`bounds_f64`.
**Fix:** borrow the shared spec and pass the tightened per-seed bounds separately (small override /
`Cow` on the bounds).

### 7. Per-eval / per-attempt allocation & recompute trims · **behaviour-preserving**
- `Problem::eval` clones the full gradient on every `gradient()` call and does a full-vector `cp == p`
  equality per cost/gradient call (`minimize.rs:91-100`). The owned-return is forced by argmin, but the
  stored clone and the equality can be trimmed (generation counter / point hash).
- `metric_gram` computes `d0[i]` (O(n²)) and discards it; `coords_from_sq` recomputes the identical
  quantity for its `< EIGVAL_TOL` reject test (`init.rs:35,64`). Return `d0` and test it directly.
- `planarity_ok` calls `improper_energy_grad(...).0`, allocating and computing a full gradient it
  discards (`checks.rs:146`). Add an energy-only improper variant.
- `stage_a_reject_energy` re-walks the chiral + 4th-dim terms already summed inside the Stage-A energy,
  and the chiral-volume loop is duplicated verbatim with `dist_geom_energy_grad` (`forcefield.rs:119-121`,
  `:79-110` vs `:206-246`). Have the DistGeom eval return the chiral+4th sub-total for reuse.

### 8. Stage-C allocates four gradient buffers per eval and merges them — **needs re-validation**
`rust/bb-embed/src/forcefield.rs:685-693`

Each of the four Stage-C terms `vec![0.0; n*3]` internally, then three are summed into the fourth — 4
heap allocs + an extra O(n) pass every Stage-C eval. Giving each term an `accum_*(…, &mut g)` form (as
the distance terms already have) is a clear win, **but** collapsing `g[c] + (gd+gi+ga)` into sequential
`+=` reassociates the FP sum, which drives the minimiser. Either keep the exact grouped summation order
while eliminating the allocations, or re-run `corpus_parity` + `stage3` to confirm the ensembles still
pass. Note the full O(n³) `symmetric_eigen` at `init.rs:77` is **not** debt: the exact solver is a ratified
correctness decision (`embed-stochastic-layer-decisions`), so a partial (Lanczos/power) solver is
out of scope by design — this is intentional, not deferred.

---

## Tier 3 — DRY / SSoT / architecture (behaviour-preserving)

### 9. `emit_inputs.rs` is a divergent copy of `build_blocks_db2` + `ligand_from_smiles`
`rust/bb-output/src/bin/emit_inputs.rs:14-89` ≡ `rust/bb-output/src/assemble.rs:68-143`

Perceive/embed, solvate + atom-count check, solv rounding, coord rounding, and the block-chunking loop
are all re-implemented in the bin; `r4` is even redefined twice in one function (`:32`, `:58`). This is
the root of finding 2. **Fix:** the bin calls a shared `build_blocks(&mol, &solv, nside) -> Vec<BlockOut{
mol2, db2 }>`; production uses the same helper and concatenates.

### 10. File-precision rounding is triplicated and decoupled from the format strings that define it
`rust/bb-output/src/assemble.rs:102-119` · `bin/emit_inputs.rs:31-45,58,79`

The authoritative precisions live in the writer format strings (`solv.rs:264,275`, `mol2_out.rs:30,76`,
`write.rs:117`); the rounding sites re-encode them as bare factors `100.0`/`10000.0`/`10.0`/`1000.0` in
two files, mapping factor↔field by hand. Any precision change to `.solv` render silently desyncs the db2
build. **Fix:** named consts (`SOLV_ATOM_DP`, `SOLV_TOTAL_DP`, `SOLV_CHARGE_DP`, `COORD_DP`, `STRAIN_DP`)
and a `round_to_file_precision(&mut TypedMol, &mut SolvFile)` on the owning types. `hierarchy::build`
knowing sibling-format decimal precision is a leaky abstraction — this removes it.

### 11. Validation layer re-implements `check.py`'s primitives instead of importing them
The `import check` pattern already works (criterion2_gate / fairness_control / tfd_recovery use it), but
a second family re-defines the same functions:
- `rmsd_matrix` (Kabsch): canonical `check.py:168`, verbatim copy `excess_validity.py:35`.
- `recovery`: canonical `check.py:116`, copy `tfd_recovery.py:22`.
- `native_spec_for`/`embed`: canonical `check.py:68`, copies in `diagnose_tfd.py:45`, `strain_compare.py:41`,
  inline in `excess_validity.py`/`bounds_strain.py`.
- **`smiles_by_id`/`mol_id` reimplemented ~6× — and the copies are subtly wrong:** the inline versions use
  `name.split("_")[0]`, which does *not* strip the `.<prot_id>` suffix that `check.mol_id` strips. They
  agree only because `seeds_100.smi` names lack the suffix — wrong on any suffixed corpus. **This is
  correctness defect C2**, not just duplication.
- The ensemble-iteration loop (`for f in ENS.glob(...): her=json.load; Hf=reshape…`) recurs near-verbatim
  in ~8 scripts.

**Fix:** a shared `validation/parity/_common.py` (or promote in `check.py`) exposing `smiles_by_id`,
`mol_id`, `CORPUS`, `native_spec_for`, `embed`, `score`, `PROD_SEED=210185`, `rmsd_matrix`, `recovery`,
`tfd_full_matrix`, `bootstrap_ci`, and an `iter_her_native(ens_dir)` generator.

### 12. Two scripts write temp spec files into the repo tree
`bounds_strain.py:43` (`validation/parity/_spec_tmp.json`) · `excess_validity.py:70` (`_ev_spec.json`,
never unlinked — leaves an untracked artifact in source). Neither is concurrency-safe — **the race is
correctness defect C4** (concurrent runs score conformers against the wrong molecule's spec).
`check.stage4()` shows the correct `tempfile.TemporaryDirectory` pattern. **Fix:** one shared `score()`
helper using `tempfile`.

### 13. Spec-diff (incl. the stereo-equivalence canonicaliser) is forked
`validation/parity/gate_spec.py:46-151` (unused by gates.sh, but cited as "the model" in `README.md:228`
/ `RIGOR.md:30`) vs `validation/hunt/differential.py:44-103` (the wired hunt gate). Two independent
implementations of the bridge-vs-native `MoleculeSpec` diff and the stereo-double-bond canonicalisation
([[stereo-representation-equivalence]]) — they can drift, so the stereo-equivalence proof could differ
between gates. **Fix:** one `spec_diff` module (`compare`, `canon_stereo`, `canon_tors`, `canon_imp`);
decide whether `gate_spec.py`'s per-molecule mode is still needed.

### 14. Scattered constants / invocations
- Container SIF path hardcoded under **two env-var names** — `BB_ORACLE_CONTAINER` (certify / corpus_writer
  / gates.sh) vs `BB_ORACLE_SIF` (`dump.sh`). The `apptainer exec --bind $REPO:$REPO $IMG` shape is
  repeated across gates.sh, `certify.py`, `corpus_writer_gate.py`.
- AMSOL paths hardcoded absolute (`corpus_writer_gate.py:20-21`) with no override — cannot run for anyone
  else or after a path move.
- `PROD_SEED = 210185` appears as a literal in ~10 scripts.
- **Fix:** a `validation/lib/env.py` with one container path (one env-var name), `apptainer_exec()`, and
  AMSOL paths with `:?must-set` guards; `PROD_SEED` in the shared module.

### 15. Rust vec3 primitives, bounds twins, magic constants
- `cross`/`dot`/`sub`/`norm3` re-implemented in `forcefield.rs`, `checks.rs`, `bounds.rs` (and `checks::sub`
  is dim-unaware while `forcefield::sub3` is dim-aware — easy to mix up). Hoist a vec3 module into `bb-core`.
- `bounds.rs` keeps near-identical f32/f64 twins (`adjust_from_coord_map[_f64]`, `coord_map_bounds[_f64]`);
  make one generic over the element type.
- Unnamed magic `0.05` = `MAX_MINIMIZED_E_PER_ATOM` (`lib.rs:121`) and the smoothing tol `0.05` written
  four times (`bounds.rs:32,36,66,74`). Promote to named consts.

---

## Tier 4 — small cleanups

(The even-length median defect that was here is a bug — see **C1** in Correctness defects, not this list.)

- `hierarchy::build` takes an unused `&SolvFile` param (`hierarchy.rs:350`) — drop it and the import.
- `clash.rs:75-102` rebuilds geometry-independent atom-type sets per conformer; precompute once.
- `mol.rs:68` `bonded_all` does an O(n) membership scan inside its BFS (fine at molecule scale; a
  `visited` bitset removes the O(n²)).
- `gen_reference.py`'s human log is parsed by regex in `speedup_report.py:22-23` — a producer format tweak
  silently yields an empty table. Have `gen_reference.py` also emit a structured per-molecule times JSON.
- `diagnose_tfd.py:89` reaches into a hardcoded `/tmp/stage3_full.log` and silently skips its headline
  correlation if absent — take the path as a CLI arg.

---

## Consolidation plan (what to create)

1. **`validation/parity/_common.py`** — corpus/id, engine (spec/embed/score/PROD_SEED/binary paths),
   geometry (rmsd/recovery/tfd/bootstrap), and the **criterion-2 gate SSoT** (`REC_TOL`, one
   `load_band()`, `failed(row, div_min)`). `criterion2_gate.py`/`fairness_control.py` shrink to drivers;
   `check.stage3()` calls the same `failed`. Absorbs findings 1, 11, 12, and the seed scatter.
2. **`spec_diff` module** shared by `gate_spec.py` and `hunt/differential.py`. Absorbs 13.
3. **`validation/lib/env.py`** — container path (one env-var name), `apptainer_exec()`, AMSOL guards.
   Absorbs 14.
4. **Rust:** shared `build_blocks(&mol,&solv,nside)` used by `emit_inputs` + `ligand_from_smiles`
   (absorbs 5, 9, 10); a `vec3` module in `bb-core` (absorbs 15); named consts.
5. **Retire/merge:** `strain_compare.py` (MMFF, self-described "poor ruler") folds into a parametrised
   `conformer_quality.py --ruler {mmff,bounds}`; `gen_reference.py` emits structured times.

---

## Resolution (2026-08-21, same day)

Every finding was executed except two, each noted with its reasoning. Behaviour was held byte-identical
where the item was tagged behaviour-preserving, and proven so by the gates + a new bit-identity unit test.

**Correctness defects:** C1 fixed (`statistics.median`, all three sites incl. the grep-found
`mutate_corpus.py`). C2/C3/C4 fixed structurally by making `check.py` the single source — `mol_id`,
the gate criterion (`REC_TOL`/`load_band`/`gate_failed`), `native_spec_for`/`embed`, `score` (tempfile),
`rmsd_matrix`, `recovery`, `PROD_SEED` — so the divergent copies no longer exist to drift.

**SSoT / DRY:** finding 11 done (`check.iter_ensembles` generator adopted by criterion2_gate / fairness_control
/ tfd_recovery; primitives single-sourced). 12 done (tempfile `score`). 13 done (`validation/lib/spec_diff.py`
— shared stereo/torsion/improper canonicalisers; both gates import them). 14 done (`validation/lib/env.py`
+ one `BB_ORACLE_CONTAINER` name everywhere; gates.sh single IMG). 9/10 done (`assemble.rs`:
`round_solv_to_file_precision` / `round_coords_to_file_precision` / `for_each_block` + named precision
consts `COORD_DP`…`STRAIN_DP`; `emit_inputs` now calls them — the gate and production share one path).
15 done (vec3 hoisted to `bb_core::vec3`; the dead f32 bounds twins **deleted** (production used only the
f64 path) rather than genericised; `MAX_MINIMIZED_E_PER_ATOM` + `SMOOTH_TOL` named).

**Hot path:** 3 done (`build_dist_pairs` once per minimization; `accum_one_pair` the single dist-term
body). 4 done (one `√d2` per violating pair). 5 done (`for_each_block` clones topology once, not per
block). 6 done (`MoleculeSpec::with_seed_bounds` shares topology, skips the two raw-bounds matrices).
7b done (`metric_gram` returns `d0`, reused by the reject test). 7c done (`improper_energy` — no gradient
alloc — via a shared `improper_geom`). 7d done (shared `chiral_vol`; the verbatim-duplicated chiral loop
is gone).

**Tier 4:** unused `&SolvFile` dropped from `hierarchy::build`; `clash` precomputes the geometry-independent
pair list once (`clash_pairs`) instead of per conformer; `bonded_all` uses a `visited` bitset; `gen_reference`
emits `_times.json` and `speedup_report` prefers it (regex fallback kept); `diagnose_tfd` takes the stage3-log
path as arg/env.

**Intentionally NOT changed, with reasoning:**
- **7a (Problem::eval clone / `cp == p`).** The memoisation is correct and the O(n) point-equality it uses is
  dominated by the O(n²) field eval it guards; the suggested generation-counter/point-hash replaces a
  provably-correct check with a scheme that can alias, for no measurable gain. The owned gradient return is
  forced by argmin's trait. Left as-is deliberately — this is a case where "optimal" favours the current code.
- **8 (Stage-C gradient buffers) — WON'T DO, and the analysis is why.** `stage_c_energy_grad` allocates 4
  gradient vectors (`torsion`→`g`, then `gd`/`gi`/`ga`) and merges `g[c] += gd[c] + gi[c] + ga[c]`, i.e.
  `g + ((gd+gi)+ga)` where each term gradient is a complete per-term sum. Eliminating the allocations is
  fundamentally at odds with byte-identity: accumulating a term directly onto a buffer that already holds a
  prior term reassociates that term's *internal* summation with the prior one (`((gd+gi₁)+gi₂+…)` instead of
  `gd+(gi₁+gi₂+…)`), which perturbs the ratified minimiser. A genuinely byte-identical rewrite must keep each
  term's complete gradient separate before combining — which is essentially the current 4-buffer code — and
  saves at most ~1 allocation while adding mutable-buffer threading that reads *worse*. The only version that
  reaches the finding's goal (≈0 allocs) is the reassociating one, and that changes output. Per the standing
  rule (everything behaviour-preserving unless a genuine behaviour fix), this is neither, so it is
  deliberately not done. A few small heap allocs per Stage-C eval are dominated by the O(pairs) math anyway.

Behaviour proof: `bb-embed` lib tests green incl. the new `pair_eval_bit_identical_to_scan` (scan vs
pre-built-pair evals byte-identical, Stage A+B, both basin regimes) and the finite-diff gradient goldens;
`corpus_parity` on 100 molecules byte-parity with live RDKit on every axis (Stage-A/B/C gradient worst
4.5e-16/3.8e-16/1.2e-13, checks 0 disagree, coord_map_bounds **exact 0.0e0**, reject energy 3.5e-15);
`dbgate` 30/30 blocks byte-identical (modulo the documented SMILES-metadata difference); `fairness_control`
still 17/94; `criterion2_gate` native 33% ≥ base2 17% → PASS; the hunt differential still 0 divergences over
132 mutants; `gate_spec` still 100/100 after the stereo-basis switch.

**Re-audit (criterion 5):** a fresh three-reviewer sweep + clippy. Clippy introduced **zero** new lints
(the ~10 remaining are the pre-existing trivial set the audit already assessed as acceptable; the one in a
rewritten file — `emit_inputs`'s `&PathBuf`→`&Path` — was fixed). The validation reviewer surfaced that
`check.py`'s OWN `stage3`/`band` still inlined the very prologue `iter_ensembles` was created to
single-source (finding 11 hadn't reached its own module) — now adopted there too, which also removed the
`210185` literal in favour of `PROD_SEED`; plus dead imports (`fairness_control` `math`/`os`, `tfd_recovery`
`json`) and an unused `part` unpack, all cleaned.

## Execution order

All items land; this is the sequence, not a filter.

1. **Correctness defects C1–C4** first — they are bugs. C1 is a one-line change; C2/C4 land with the
   validation consolidation and C3 with it too, but the fixes are not contingent on the refactor's scope.
2. **Hot-path 3–6** (behaviour-preserving) — the embed-loop wins; re-run `gates.sh fast` + `corpus_parity`
   after the FF-loop change to confirm byte-identity.
3. **Validation consolidation** (findings 11, 12, 14) — the shared `_common.py` / `spec_diff` / `env.py`
   modules; single-sources the criterion-2 gate, the geometry primitives, and the container/AMSOL/seed
   constants. This is where C2/C4 land and the duplication behind C3 is removed.
4. **Output-pipeline single-sourcing** (findings 9, 10) — `build_blocks` shared by `emit_inputs` +
   `ligand_from_smiles`; file-precision constants and rounding onto the owning types.
5. **Rust DRY / constants** (finding 15) — the `bb-core` vec3 module, the f32/f64 bounds twins made
   generic, the named consts.
6. **Finding 8** (Stage-C gradient buffers) — done, gated by a full parity re-validation since it
   perturbs the numerics; keep the exact summation order if the reassociation doesn't pass.
