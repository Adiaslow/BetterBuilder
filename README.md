# BetterBuilder

Macrocycle ligand building for DOCK, in Rust. A reimplementation of the `build_ligands.py` pipeline
(the oracle): conformers by the RDKit ETKDGv3 core-pin recipe, AMSOL solvation, torsion strain, and
the mol2 + db2 tarball DOCK consumes — intended as a drop-in replacement with the same output
behaviour, built from different machinery. The production pipeline uses no RDKit.

## Known differences from `build_ligands.py`

- **Full conformer count.** The oracle silently returns fewer conformers when no sidechain can be
  embedded around a core; BetterBuilder replaces that core, so every molecule gets the recipe's full
  count.
- **Atom order.** The oracle rewrites each SMILES in RDKit canonical form before building, so its files
  number atoms in canonical order; BetterBuilder keeps the input SMILES' order. Same molecules; atom
  indices differ between the two systems' files.
- **Protomer ids and order.** Input lines may give a protomer id (`smiles name prot_id`). A repeated id
  takes the next free id above it (the oracle's renumbering can still collide on three repeats), and
  molecules keep input order (the oracle sorts them by `name.prot_id` text, which also orders its
  archive).
- **Archive.** The tarball is finalized properly (the oracle closes its archive only after moving it,
  which can leave the moved copy truncated) and has no `.dock_version` member.
- **Names.** ZINC-22 ids carrying a tranche are filed in their tranche directory, per the ZINC-22
  convention; legacy `ZINC00…` ids go to the top level (the oracle files them under an `H00` tranche
  that does not exist). Any molecule whose name cannot name a file (empty, `.`, `..`, or containing
  `/`) is reported and skipped.

## Architecture

| Crate (`rust/`) | Role |
|---|---|
| `bb-perceive` | Molecular perception reproducing RDKit's: SMILES parsing, rings, aromaticity, valence, AddHs, stereo, SMARTS matching, experimental torsions. |
| `bb-bounds` | The distance-geometry bounds matrix (RDKit's `setTopolBounds`) and the force-field parameters it reads. |
| `bb-spec` | Setup: assembles a molecule's `MoleculeSpec` from `bb-perceive` and `bb-bounds`. CLIs `bb-spec-native`, `bb-spec-native-batch`. |
| `bb-core` | The `MoleculeSpec` data contract shared by setup and the engine, and the triangle smoothing both use. |
| `bb-embed` | The engine: staged distance-geometry minimization (Stages A/B/C), acceptance checks, and the two-stage core-pin recipe. CLIs `bb-embed`, `bb-batch`. |
| `bb-solv` | Solvation: conformer → MOPAC Z-matrix → AMSOL 7.1 inputs, invocation, and the `.solv` output. AMSOL itself is external and configured. |
| `bb-strain` | Torsion-strain energies (a port of the pipeline's `Torsion_Strain`), carried in the db2. |
| `bb-db2` | The DOCK `.db2` format: data model, parser, byte-exact writer. |
| `bb-output` | The output artifacts — typed molecule, mol2, db2 hierarchy — and the tarball. CLI `bb-build` runs the whole pipeline. |
| `bb-py` | PyO3 bindings: the `betterbuilder` Python module (setup + embed). |
| `bb-rdkit` | The bridge to the patched RDKit C++ (`cxx` FFI). Not in the production pipeline: it serves the parity gates, fixture generation, and the bridge spec CLI `bb-spec`. |

The archived GPU implementation and one-off prototypes live in `.archive/` (gitignored).

## Build

The production tools are pure Rust and build with cargo alone:

```sh
cd rust && cargo build --release -p bb-output -p bb-spec -p bb-embed   # bb-build, bb-spec-native, bb-embed
```

The full workspace (`cargo build --release` in `rust/`) also builds `bb-rdkit`, which links a
**patched RDKit build** (2026.09.1pre + the amide-torsion patch) that is not bundled. To build that
dependency and everything it needs into `$HOME/toolchain`:

```sh
toolchain/build-all.sh                         # GCC, Eigen, Boost, patched RDKit — no root needed
. toolchain/env.sh                             # exports the prefixes the build reads
cd rust && cargo build --release
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

The whole pipeline, SMILES file → DOCK tarball:

```sh
bb-build input.smi output.tar.gz workdir [seed]
# input.smi: one `smiles name [prot_id]` per line; AMSOL from BB_AMSOL_EXE (+ BB_AMSOL_LD_LIBRARY_PATH,
# BB_AMSOL_TIMEOUT, BB_AMSOL_MAX_CONCURRENT, or a BB_CONFIG file). A failing line or molecule is
# reported on stderr and skipped.
```

Conformers only (setup → embed):

```sh
bb-spec-native "C1CCOCC1" | bb-embed - recipe 210185 > conformers.json
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

`rust/gates.sh` is the gate suite's entry point: `fast` (toolchain-free gates against committed
goldens), `verify` (regenerate the goldens from live RDKit and byte-diff them), `certify` (our RDKit
build against her container's), `live` (RDKit-parity gates on a chosen corpus), `hunt` (adversarial
macrocycles, native vs bridge), `dbgate` (db2 writer vs her `mol2db2`, byte-identical) and `bench`
(each Rust component against RDKit doing equal work). `validation/parity/` compares BetterBuilder
against the oracle stage by stage — see [`validation/parity/README.md`](validation/parity/README.md).

## Repository layout

```
rust/              Rust workspace — the crates above, and gates.sh
validation/        parity harness, gate fixtures and corpora (validation/README.md)
rdkit-patch/       the patch defining the patched-RDKit build dependency (BB_RDKIT_ROOT)
toolchain/         scripts building that dependency + its compiler and libraries into $HOME/toolchain
docs/              the 2026-08-21 code audit
```

External reference kept on disk but **out of the repo** (gitignored — not the BetterBuilder deliverable):
`oracle/` + `oracle_local/` (the RDKit oracle — validation/benchmark source of truth), `old_code/`
(the DOCK downstream pipeline the output stage reproduces), `.archive/` (GPU impl + prototypes).

## Status

`bb-build` implements the full pipeline. Stage-by-stage parity with the oracle is recorded in
[`validation/parity/`](validation/parity/README.md); measurements made before the 2026-10 changes are
marked there for re-measurement. Not yet done:

- **End-to-end gate** — compare the tarball against the oracle's on her container, with AMSOL.
- **Throughput** — previously reported speed figures (~3.4–4.5× per core) predate the current code;
  re-measure against the oracle.

[maturin]: https://www.maturin.rs
