# BetterBuilder

Fast, faithful macrocycle conformer generation. A pure-Rust reimplementation of the RDKit ETKDGv3
core-pin recipe used by the `build_ligands.py` oracle — validated to reproduce the oracle's conformer
ensembles at the run-to-run noise floor, while running several× faster per core.

## Architecture

Perception (**setup**) and minimization (**embed**) are split cleanly; only the embed — the ~91% of
the work — is reimplemented, and setup is borrowed from RDKit:

| Crate (`rust/`) | Role |
|---|---|
| `bb-core`  | The `MoleculeSpec` data contract (bounds matrix, torsions, chiral sets, impropers) shared across crates. |
| `bb-rdkit` | **Setup**: SMILES → `MoleculeSpec` via `cxx` FFI to the patched RDKit C++. Panic-free library; the only piece borrowed from RDKit. |
| `bb-embed` | **The engine**: staged distance-geometry minimization (Stages A/B/C), acceptance checks, the two-stage core-pin recipe. Pure Rust, best-in-class numerics (argmin, nalgebra). |
| `bb-py`    | Thin [PyO3] bindings → the `betterbuilder` Python module. |

The two CLIs pipe setup into embed over the `MoleculeSpec` JSON contract. The archived GPU
implementation and one-off prototypes live in `.archive/` (gitignored).

## Build

BetterBuilder links a **patched RDKit build** (2026.09.1pre + the amide-torsion patch) that is *not*
bundled — you provide its location via the environment (the build fails with a clear message if unset).
See [`rdkit-patch/`](rdkit-patch/README.md) for the patch and how to build that dependency.

```sh
export BB_RDKIT_ROOT=/path/to/patched/rdkit    # the dir containing Code/ and lib/  (RDKit's RDBASE also works)
# BB_DEPS_PREFIX auto-probes /opt/homebrew, /usr/local, /usr for boost + eigen — override if needed.

cd rust && cargo build --release               # engine + CLIs (bb-spec, bb-embed, bb-batch)
```

The Python module (needs [maturin]):

```sh
cd rust/bb-py && maturin develop --release      # builds + installs the `betterbuilder` wheel
```

`rust/.cargo/config.toml` sets `target-cpu=native` for the build machine (guard-validated
fidelity-neutral). For a binary shipped to different hardware, override with the deploy arch, e.g.
`RUSTFLAGS="-C target-cpu=x86-64-v3"`.

## Use

Command line (setup → embed):

```sh
bb-spec "C1CCOCC1" | bb-embed - recipe 210185 > conformers.json
# {"n_atoms": N, "conformers": [[x0,y0,z0, x1,…], …]}   (n_atoms*3 per conformer, RDKit atom order)
```

Python:

```python
import betterbuilder
result = betterbuilder.embed("C1CCOCC1", seed=210185)
result.n_atoms          # heavy atoms + H, RDKit (SmilesToMol + AddHs) order
result.conformers[0]    # flat [x0, y0, z0, x1, …] for conformer 0
len(result)             # number of conformers (chosen by the recipe from the rotatable-bond count)
```

## Validate

```sh
<python-with-rdkit> validation/regression_guard.py     # exit 0 = PASS, 1 = REGRESSION
```

Runs the candidate (fixed seed → deterministic) on 100 reference molecules, compares its conformer
ensembles to the oracle's by heavy-atom nearest-neighbour RMSD, and gates against a pinned baseline
with noise-derived thresholds. Current: **median NN 1.260 Å**, sitting at the oracle's own
seed-to-seed noise floor. See `validation/README.md`.

## Repository layout

```
rust/              Rust workspace — the engine (bb-core, bb-rdkit, bb-embed) + Python bindings (bb-py)
validation/        regression guard, pinned baseline, oracle reference ensembles
rdkit-patch/       the patch defining the patched-RDKit build dependency (BB_RDKIT_ROOT)
docs/              architecture, validity methodology, status
```

External reference kept on disk but **out of the repo** (gitignored — not the BetterBuilder deliverable):
`oracle/` + `oracle_local/` (the RDKit oracle — validation/benchmark source of truth), `old_code/`
(DOCK downstream-pipeline reference for the output-gap work), `.archive/` (GPU impl + prototypes).

## Status

Scientifically complete: fidelity validated at the ensemble noise floor, ~3.4–4.5× faster per core
than the oracle (widening with molecule size), deterministic, `clippy`-clean, panic-free setup library.
Remaining before production deployment — see `docs/03-status.md`:

- **Output gap** — emit the SDF/db2 the DOCK pipeline consumes (the `betterbuilder` Python API is the vehicle).
- **Full-set validation** — the guard samples 100 of 1000 molecules; the rest need oracle references generated.
- **Linux/x86 build** of the patched RDKit + a deploy-arch guard re-run (FMA rounding).

[PyO3]: https://pyo3.rs
[maturin]: https://www.maturin.rs
