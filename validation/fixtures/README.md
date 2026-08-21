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

## `divya/` — embed-ensemble + output goldens (her seeded pipeline) — TODO

The authoritative oracle. To be generated from `build_ligands.py` with her seeds (`randomSeed =
210185 + j`) via her container, and compared **distributionally** (RMSD/TFD) because native's
RNG/eigensolver/minimizer differ from hers, so identical seeds do not give identical conformers.
Note: her solvation embed (`build_ligands.py` ~line 191) appears **unseeded**, so its charge branch
is reproducible only as a one-sample capture, not "deterministic given seeds" — verify per-stage
seeding when generating this tier. (The previous `validation/reference/*.sdf` were stale and removed.)
