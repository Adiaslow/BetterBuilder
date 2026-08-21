# Validation rigor rubric

The bar every gate must clear before it counts as validated. Written down so rigor is *enforced*,
not rediscovered per gate. A gate that does not meet its row here is provisional, and must say so.

## Every gate, regardless of type

- **Names and characterizes its corpus.** A gate reports only on the chemistry it ran against. The
  corpus is an explicit input with a `.provenance` file, and is characterized against the property
  under test before use — for this pipeline that is Divya's amide reweighting (fires on 83% of the
  production workload; the only reason her RDKit is not stock). A gate measured on a corpus where
  the patch fires at 16% has not been measured on the workload. (See `README.md` "Choosing a corpus".)
- **Selection is not done by our own code where the oracle can do it.** Picking candidates with our
  own code samples only where it already agrees with itself; a molecule it wrongly fails to nominate
  never gets looked at.
- **Reference values come only from the oracle executing** — never BetterBuilder, never a modified
  copy of the oracle. (See `README.md` "The oracle".)

## Deterministic gates (spec, bounds, FF gradients, checks, db2, solvation-on-fixed-geometry)

The output is a function of the input; the oracle gives one right answer.

- **Exact or derived-tolerance match.** Byte-identical where the format is exact; otherwise a
  tolerance *derived* from a named cause (f32 storage → 2e-6 Å; machine precision → ~1e-12), never a
  round number picked to pass.
- **Run over the provenanced corpus, not a hand-picked set.** A formula/logic gate on a curated
  handful proves the formula *on that handful*. Coverage still bites: the `coord_map_bounds` bug
  showed only on a two-ring molecule; benzoic passed. Corpus scale is not optional just because the
  comparison is exact.
- **Report worst-case magnitude and a disagreement count**, not a boolean. `gate_spec` is the model:
  every field, every molecule, worst |Δ|, count of molecules differing.
- **Cover the geometry, not just the molecule.** For gates evaluated at coordinates (FF, checks),
  include the degenerate / near-collinear cases as well as embedded conformers — those are what
  caught the dihedral and angle-constraint divergences.

## Stochastic gates (conformer ensembles; eventually solvation-from-our-conformer)

The embed's RNG/minimizer differ from hers by construction, so the output cannot byte-match; only
distributions can be compared.

- **Distance is symmetry-corrected.** RMSD minimized over the heavy-atom automorphism group (a
  symmetric substituent must not be penalized for relabelling), *and* TFD (torsion-fingerprint
  deviation) — the ring-aware, symmetry-aware convention for macrocycles, where Cartesian RMSD
  washes out ring/rotamer differences. Report both; where they disagree is diagnostic.
- **Baseline is the oracle's own run-to-run spread, not ours.** An implementation that samples too
  narrowly shows a small self-spread and scores well against itself. Normalize to *her* spread. The
  gold standard is her spread across **independent seeds** (run her twice); her within-run
  seed-to-seed spread (first half vs second half of one ensemble) is the interim baseline.
- **Two-sided.** Nearest-neighbour coverage alone rewards collapse into one basin. Guard diversity
  too — within-seed (sidechain sampling around a pinned core) and between-seed (core conformations).
  Matching her spread is the requirement; being more homogeneous is a defect, not an achievement.
- **Per-molecule AND aggregate, with uncertainty.** Report the per-molecule distribution and a
  per-molecule pass count, and a bootstrap CI on the aggregate. Gate on the CI, not a point
  estimate — an aggregate within band while 26% of molecules fail is not a pass.
- **Thresholds are derived, not picked.** The acceptance band comes from the oracle's own
  independent-seed variability, measured. A ±15% placeholder is provisional and must be labelled so.
- **Report an interpretable coverage curve** (fraction of her conformers recovered within τ) beside
  the normalized ratio, and always beside *her own* self-recovery.
- **Include an energetic axis.** Geometry coverage without energetics is half the picture; compare
  the strain (or a consistent-FF energy) distribution, not only shape.

## Running the suite

One committed entry point (`rust/gates.sh`), sourcing the committed `toolchain/env.sh` (never a
scratchpad copy), with three tiers:

```
rust/gates.sh fast     # native correctness + parity vs committed goldens — NO RDKit toolchain
rust/gates.sh verify   # regenerate goldens from live RDKit, byte-diff committed (anti-staleness)
rust/gates.sh live     # live RDKit-parity gates + the corpus-scale gate (seeds_5000)
rust/gates.sh          # all three
```

A red gate aborts non-zero: red is a signal to fix the code, never licence to weaken or delete a gate.

### Fixture tiers (`validation/fixtures/`, see its README)

Parity goldens are a **cache of an oracle executing** — generated only by the live bridge, never
hand-authored — with each record's input stored beside its output so the byte-identical-input
coupling (a diff = a formula difference) is preserved.

- `rdkit/` — deterministic-component goldens = Divya's **patched** RDKit executing (her pipeline
  delegates these pieces to it). Committed on `seeds_100`; `gen-fixtures` regenerates them; the
  `fast` tier validates native against them with no toolchain, and `verify` byte-checks them against a
  fresh live run (pinned to `VERSION`). The full `seeds_5000` authority stays the live `corpus_parity`.
- `divya/` — embed-ensemble + output goldens = her **seeded** pipeline, compared distributionally
  (RMSD/TFD). Authoritative, TODO (generate from `build_ligands.py`).

Gates evaluated at a **native-generated** geometry (Stage-C at native's conformer, coord_map at
native's pins) cannot be pre-serialized — their input is native's own output — so they stay live.

## Standing status (fill in as gates are brought to the line)

| Gate | Type | Corpus | To standard? |
|---|---|---|---|
| `gate_spec` (all spec fields) | deterministic | 5000 seeds | yes |
| metric-matrix embedding (`bb-embed/tests/embedding.rs`) | native-correctness | random point clouds (ground truth) | yes — recovers embeddable geometry to 1.3e-12 Å |
| Stage-A reject energy (`corpus_parity` reject-energy axis) | deterministic (RDKit parity) | 5000 seeds | yes — `stage_a_reject_energy` == RDKit `calcEnergy` to ~1e-12 rel over the corpus (the 0.05/atom reject decision is identical by construction when the energy matches) |
| bounds 12/13/14/15, raw_bounds; rings; valence; AddHs; hybridization; angles | deterministic | 600–5000 | yes |
| FF gradients (Stage A/B/C) | deterministic | 5000 (`corpus_parity`, 83% patch) | yes — A/B 6e-16, C 7e-7 (proven V=100-torsion 1/sinφ FP floor, tol 1e-6) |
| embed checks (6) | deterministic | 5000 | yes — 0/5000 disagree |
| `coord_map_bounds` | deterministic | 5000 | yes — bit-identical to RDKit (0.0) |
| Stage 4 solvation (fixed geometry) | deterministic | 24 | partial — coverage; 4b unmeasured |
| Stage 3 ensembles | stochastic | 23 | metric hardened; **thresholds provisional (Phase 2 step 4)** |
| strain / energy | — | — | **not built** |
| aromaticity / kekulization | assumption | corpus-scoped | gated only as-written |
