# Oracle-provenance certification

**CONFORMANCE.md Phase 1.** Every deterministic golden under `validation/fixtures/rdkit/` is only
Divya-faithful if "our from-source RDKit build == her actual production build." This tier *proves*
that, rather than trusting our patch reconstruction (`rdkit-patch/`).

## Method

Divya's production container (`build_macrocycle_final.sif`) ships her real RDKit **libs**
(`/soft/rdkit_libs`, `2026.09.1pre`, carrying her amide patch) and a C++ compiler. Our host toolchain
(RDKit **headers** of the same version, boost 1.82.0 = hers, eigen) is bind-mounted into the
container. So we compile a small **oracle dumper** — a standalone mirror of the corresponding
`bb-rdkit/src/bridge.cc` function — with our headers against **her libs**, run it in-container, and
diff its output against our host-bridge golden at the *identical* stored geometry. Same RDKit version
+ same boost → ABI-compatible; a difference is a build/patch difference between her RDKit and ours.

Result: **byte-identical (`0.00e+00`) over the 100-molecule corpus for all six deterministic oracles**
(Phase 1 complete) — the Stage-C match is the amide-patch proof (her patched 80.0/130.9 torsion
constants reproduced exactly).

## Dumpers (mirrors of `bridge.cc`; that file is the single source of truth)

- `stage_a_dump.cc` — `stage_a_ff_grad` (`constructForceField`: distance + chiral + 4th-dim). Weights
  are argv, so it certifies **both Stage-A** (firstMinimization, 1.0/0.1) **and Stage-B**
  (minimizeFourthDimension, 0.2/1.0) — same FF code, different weights.
- `stage_c_dump.cc` — `stage_c_ff_grad` (getExperimentalTorsions incl. the amide patch + impropers +
  constraints). **The amide-patch certification.**
- `checks_dump.cc` — `embed_checks`: the six post-embed accept/reject checks (EmbeddingOps + the
  planarity FF), with the `detail::EmbedArgs` struct forward-declared byte-for-byte.
- `coord_map_dump.cc` — `coord_map_bounds` (`adjustBoundsMatFromCoordMap` + `triangleSmoothBounds`).
  Certifies **coord_map** with embedded-conformer pins and the plain smoothed **bounds** at k=0 pins.

## Running

```
rust/gates.sh certify        # compiles the dumpers in-container against her libs, diffs vs goldens
```

`gates.sh certify` compiles each dumper with the pinned include paths (see `bb-rdkit/build.rs`) and
runs `certify.py <N> <dumper> <golden.jsonl>`. Env `BB_ORACLE_CONTAINER` overrides the `.sif` path.
A non-byte-identical result fails the tier — it means our from-source build has drifted from hers
(regenerate/rebuild), which is a signal to fix, never to loosen.
