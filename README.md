# BetterBuilder

Fast macrocycle conformer generation. A pure-Rust reimplementation of the RDKit ETKDGv3 core-pin
recipe used by the `build_ligands.py` oracle, intended as a drop-in replacement: the same output
behaviour, different and faster machinery.

Parity with the oracle is verified stage by stage for perception, recipe parameters and conformer
ensembles; solvation is not yet compared, and strain/db2 are not built, so the end-to-end tarball
comparison does not exist yet. See [`validation/parity/`](validation/parity/README.md).

## Architecture

Perception (**setup**) and minimization (**embed**) are split cleanly; only the embed — the ~91% of
the work — is reimplemented, and setup is borrowed from RDKit:

| Crate (`rust/`) | Role |
|---|---|
| `bb-core`  | The `MoleculeSpec` data contract (bounds matrix, torsions, chiral sets, impropers) shared across crates. |
| `bb-rdkit` | **Setup**: SMILES → `MoleculeSpec` via `cxx` FFI to the patched RDKit C++. Panic-free library; the only piece borrowed from RDKit. |
| `bb-embed` | **The engine**: staged distance-geometry minimization (Stages A/B/C), acceptance checks, the two-stage core-pin recipe. Pure Rust, best-in-class numerics (argmin, nalgebra). |
| `bb-solv`  | **Solvation**: conformer coordinates → MOPAC Z-matrix → AMSOL7.1 input files, invocation, and the `.solv` output the DOCK pipeline consumes. AMSOL itself is external and configured. |
| `bb-py`    | Thin [PyO3] bindings → the `betterbuilder` Python module. |

The two CLIs pipe setup into embed over the `MoleculeSpec` JSON contract. The archived GPU
implementation and one-off prototypes live in `.archive/` (gitignored).

## Build

`bb-core`, `bb-embed` and `bb-solv` are pure Rust and build with cargo alone. `bb-rdkit` links a
**patched RDKit build** (2026.09.1pre + the amide-torsion patch), which is not bundled.

To build that dependency and everything it needs into `$HOME/toolchain`:

```sh
toolchain/build-all.sh                         # GCC, Eigen, Boost, patched RDKit — no root needed
. toolchain/env.sh                             # exports the prefixes the build reads
cd rust && cargo build --release               # engine + CLIs (bb-spec, bb-embed, bb-batch)
```

See [`toolchain/`](toolchain/README.md) for prerequisites and per-package options, and
[`rdkit-patch/`](rdkit-patch/README.md) for what the patch changes. To use an RDKit built elsewhere,
set `BB_RDKIT_ROOT` (the directory containing `Code/` and `lib/`; `RDBASE` also works) plus
`BB_BOOST_PREFIX` and `BB_EIGEN_PREFIX`, or `BB_DEPS_PREFIX` when boost and eigen share one prefix.
`BB_CXX_PREFIX` names the C++ runtime RDKit was built against, when it is not the system one.

The Python module (needs [maturin]):

```sh
cd rust/bb-py && maturin develop --release      # builds + installs the `betterbuilder` wheel
```

`rust/.cargo/config.toml` sets `target-cpu=native` for the build machine. For a binary shipped to
different hardware, override with the deploy arch, e.g.
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

`validation/parity/` compares BetterBuilder against the oracle's unmodified code, stage by stage —
perception, recipe parameters, conformer ensembles — using values read from her code executing in her
production container. `regression_guard.py` is the older ensemble-level gate; it has no reference
data at present. See [`validation/parity/README.md`](validation/parity/README.md).

## Repository layout

```
rust/              Rust workspace — the engine (bb-core, bb-rdkit, bb-embed) + Python bindings (bb-py)
validation/        parity harness (validation/parity/) + molecule lists
rdkit-patch/       the patch defining the patched-RDKit build dependency (BB_RDKIT_ROOT)
toolchain/         scripts building that dependency + its compiler and libraries into $HOME/toolchain
docs/              architecture, validity methodology, status
```

External reference kept on disk but **out of the repo** (gitignored — not the BetterBuilder deliverable):
`oracle/` + `oracle_local/` (the RDKit oracle — validation/benchmark source of truth), `old_code/`
(DOCK downstream-pipeline reference for the output-gap work), `.archive/` (GPU impl + prototypes).

## Status

Engine complete and deterministic; `clippy`-clean, panic-free setup library. Parity with the oracle
is verified for perception, recipe parameters and conformer ensembles
([`validation/parity/`](validation/parity/README.md)). Previously reported speed figures
(~3.4–4.5× per core) predate the current setup and need re-measuring.
Remaining before production deployment — see `docs/03-status.md`:

- **Output gap** — strain and db2, so the pipeline can emit the tarball DOCK consumes.
- **End-to-end gate** — compare that tarball against the oracle's, which needs the above.
- **Throughput** — re-measure against the oracle once the pipeline is complete.

[PyO3]: https://pyo3.rs
[maturin]: https://www.maturin.rs
