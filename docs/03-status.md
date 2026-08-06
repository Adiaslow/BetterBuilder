# BetterBuilder — status

Snapshot of where the work stands. Architecture in `docs/00`; validity background (parked) in `docs/01`.

## Done + validated

| Piece | State |
|---|---|
| **Setup** (`bb-rdkit`): SMILES → `MoleculeSpec` via patched-RDKit FFI | Panic-free `Result` API; carries the patched amide torsions (proven via FFI). ✅ |
| **Engine** (`bb-embed`): staged embed (A/B/C), acceptance checks, core-pin recipe | Every step matched to RDKit source; finite-diff gradient checks; triangle-smooth bit-exact vs RDKit. ✅ |
| **Fidelity** (`validation/regression_guard.py`) | Median NN **1.260 Å** over 100 molecules — at the oracle's seed-to-seed noise floor. Deterministic. ✅ |
| **Speed** (vs the real `build_ligands.local.py` oracle) | **~3.4–4.5× per core**, widening with molecule size; + `lto`/`target-cpu` build flags (guard-neutral). ✅ |
| **Python API** (`bb-py`): `betterbuilder.embed(smiles, seed)` | Thin PyO3; **bit-for-bit matches the CLI engine**; invalid SMILES → `ValueError`. ✅ |
| Hygiene | `clippy`-clean, rustfmt, workspace metadata, env-configured build (no hardcoded paths). ✅ |

## How to run (from `rust/`, with `BB_RDKIT_ROOT` set — see the root README)

```bash
cargo build --release                                   # engine + CLIs
cargo test --release                                    # unit + finite-diff + discriminating checks tests
(cd bb-py && maturin develop --release)                 # the betterbuilder Python module
../validation/regression_guard.py                       # fidelity guard (needs a python with RDKit)
```

## Remaining before production deployment

1. **Output gap** — emit the SDF/db2 the DOCK pipeline consumes; the `betterbuilder` Python API is the
   vehicle, `old_code/strain/` the downstream reference.
2. **Full-set validation** — the guard samples 100 of 1000 molecules; generate oracle reference
   ensembles for the rest (hours of oracle compute) and re-run.
3. **Linux/x86** — build the patched RDKit on the deploy architecture; re-run the guard there (FMA
   rounding differs from the M2/ARM dev build).
4. **`checks.rs`** — 3 discriminating unit tests added (chiral/tetrahedral/stereo rejection); planarity
   is still only guard-covered at the ensemble level.
5. **Scale** — parallel efficiency validated to 12 cores (7.4×); the target is 500 (NUMA/memory-bandwidth
   behavior untested at that scale).

## Archived (`.archive/`, gitignored; details in project memory)

- **GPU pipeline** (CubeCL→Metal): correct and guard-passing, but retry-bound and *loses* to the CPU
  engine end-to-end; residual rare intra-kernel race needs CUDA `racecheck` to pin down.
- **Rejected CPU accelerations** (prototyped and measured): SoA-over-conformers SIMD (~1.7× but requires
  FIRE, whose retry blowup dominates), spatial active-set (~40% of active constraints are far upper-bound
  pairs a spatial hash can't find), f32 force field (~1.05× on the large molecules that matter).
