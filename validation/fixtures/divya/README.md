# Divya reference ensembles — the Phase-3 distributional oracle

The authoritative oracle for conformance target **B** (CONFORMANCE.md): native's conformer *ensemble*
is compared **distributionally** to Divya's actual seeded ensembles (never conformer-for-conformer,
since native's RNG/eigensolver/minimizer differ from hers by design).

## `gen_reference.py` — her ensemble, from her container

A **verbatim** replica of the ensemble block in her `build_ligands.py`: the seeded two-stage core-pin
recipe —
1. seeded core embed (`ETKDGv3`, `useMacrocycleTorsions`, `useExpTorsionAnglePrefs`, `randomSeed =
   210185 + j`) for `rdkit_confs_1` core conformers;
2. pin the largest ring + exocyclic `=O` + exocyclic C on ring N (`cmap`);
3. seeded sidechain expansion (`EmbedMultipleConfs`, `useRandomCoords`, `randomSeed = 210185`,
   `SetCoordMap(cmap)`) — `sidechain_confs` per core.

`(rdkit_confs_1, sidechain_confs)` is `(20, 10)` if `count_exo_rotatable ≤ 2` else `(10, 20)` → ~200
conformers. Only the surrounding pipeline and I/O are replaced (corpus-read + SDF-write); the embed
logic is copied, not reimagined — fidelity to her output is the point, so do not "improve" it.

Run inside her production container so it uses **her actual RDKit build** (2026.09.1pre + her amide
patch — the same one Phase-1 certifies our build byte-equal to):

```
apptainer exec --bind <repo>:<repo> <build_macrocycle_final.sif> \
  python3 validation/fixtures/divya/gen_reference.py <corpus.smi> validation/fixtures/divya/ensembles
```

Writes `ensembles/<name>.sdf` (all conformers) per molecule.

## Why not committed

Verified: each ensemble is ~1.4 MB (~138 MB for `seeds_100`) and **byte-deterministic across runs**
(all seeds fixed). So the ensembles are **regenerated on demand** into the gitignored `ensembles/`
dir — never committed — exactly as the certify tier recompiles dumpers rather than committing
binaries. Determinism means a regenerated set is identical to any prior one.

## Acceptance band (the distributional gate, next)

The band is her *own* variability, not a picked tolerance. Her single-seed-base ensemble is
deterministic, so the between-independent-seed spread is measured by regenerating at a **second seed
base** (e.g. `210185` vs another) and comparing; the within-ensemble core-seed-to-core-seed spread is
the interim baseline. Native (`bb_embed::embed_recipe`) is then required to sit within that band on
symmetry-corrected RMSD + TFD + strain, per molecule + bootstrap CIs (RIGOR.md stochastic gate).
