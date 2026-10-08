# Parity fixtures

Serialized oracle outputs that native gates validate against, so the fast inner loop needs **no live
oracle**. Fixtures are a *cache of an oracle executing*, never hand-authored, and each record carries
its input alongside the output so the byte-identical-input coupling is preserved (a diff means a
formula difference, not a plumbing one).

## `rdkit/` — deterministic-component goldens (her patched RDKit)

Divya's pipeline delegates the deterministic pieces (bounds, force field, checks, coord-map) to the
same patched RDKit our bridge runs, so "native == RDKit here" *is* "native == what her pipeline
computes here". These goldens exist to make the component gates fast and toolchain-free, and to
**localize** any end-to-end divergence to a stage.

- `VERSION` — the patched-RDKit version the goldens were generated against (bump-guard).
- `stage_a_grad.jsonl` — `{smiles, coords[n*4], grad[n*4]}`: RDKit's Stage-A force-field gradient at a
  fixed 4D geometry. (More gate files added as they are converted.)

Generated **only** by the live bridge:

```
cargo run -p bb-rdkit --bin gen-fixtures -- validation/seeds_100.smi validation/fixtures/rdkit
```

Consumed by `bb-embed/tests/parity_goldens.rs` (no RDKit toolchain). The committed goldens are on
`seeds_100`; the full `seeds_5000` authority stays a **live** `gates.sh --verify` pass that
regenerates and byte-diffs these files, so they cannot silently rot and a version bump is caught.

## `zinc22/` — ZINC-22 id → 3D-database directory (the ZINC maintainers' code)

`tranche_dirs.tsv` — `tranche, sub_id, zinc_id, directory`: ZINC-22 ids produced by the maintainers' own
encoder and the directory their own `get_zinc_directory_hash` files each under. Covers every tranche in
their list of real 3D tranches plus every logP bin at heavy-atom counts 0, 1, 29, 30 and 61. Generated
**only** by running their code at pinned commits (the URLs are in the file's header):

```
python3 validation/fixtures/zinc22/gen_tranche_dirs.py > validation/fixtures/zinc22/tranche_dirs.tsv
```

Consumed by `bb-output`'s `zinc22` tests. Legacy ids (`ZINC00…`, ZINC-20 and earlier) are not in it:
the ZINC-22 numbering defines them as carrying no tranche, which their directory function does not
honour, so those cases are tested against that definition instead.

## `divya/` — her seeded conformer ensembles

`divya/gen_reference.py` regenerates her seeded ensembles (`randomSeed = 210185 + j`) with a verbatim
replica of the ensemble block of `build_ligands.py`; see [`divya/README.md`](divya/README.md). The
ensembles are byte-deterministic and large, so they are regenerated on demand into the gitignored
`divya/ensembles/` and never committed. They are compared **distributionally** (RMSD/TFD), because
native's RNG/eigensolver/minimizer differ from hers, so identical seeds do not give identical
conformers. Her solvation embed (`build_ligands.py:191`) is **unseeded**, so its charge branch is
reproducible only as a one-sample capture, not "deterministic given seeds". (The previous
`validation/reference/*.sdf` were stale and removed.)
