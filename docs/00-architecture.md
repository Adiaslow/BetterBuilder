# BetterBuilder — architecture

A pure-Rust reimplementation of the macrocycle ligand builder's **conformer-generation embed** — the
ETKDG stage that is ~91% of build time. It reproduces RDKit's *algorithm* faithfully (validated by
ensemble RMSD to the RDKit oracle at the run-to-run noise floor — not bit-exact) while running several×
faster per core. Setup is borrowed from RDKit; only the embed is reimplemented.

## The shape: RDKit setup → Rust engine → conformers

Perception (**setup**) is cheap and not GPU/CPU-shaped; the embed is expensive and is the ~91%. So they
are split across a numeric contract:

```
┌ setup (bb-rdkit, RDKit FFI, once) ─────────┐   ┌ embed engine (bb-embed, per conformer) ─────────┐
│ SMILES → perception → bounds matrix +      │ → │ init → Stage A (4D) → Stage B (squeeze) →       │
│ exp-torsions + chiral/improper/constraints │   │ Stage C (3D torsions) → acceptance checks/retry │
│   ⇒ MoleculeSpec (numeric arrays)          │   │   ⇒ conformer coordinates                        │
└────────────────────────────────────────────┘   └──────────────────────────────────────────────────┘
```

**Perception boundary (load-bearing):** RDKit perceives once and hands the engine a numeric
`MoleculeSpec`; the engine does only per-conformer geometry and never perceives. RDKit's setup is
exposed and cheap, so reimplementing it (re-deriving SMARTS / torsion libraries) would be the
anti-goal — see `docs/01` and the oracle methodology.

`MoleculeSpec` (the setup→embed contract, `bb-core`):

```
n_atoms, dim (3/4)
bounds:        N×N f32  (row<col = upper bound, row>col = lower)   — getMoleculeBoundsMatrix
exp_torsions:  [{ atoms[4]; V[6]; signs[6] }]                       — getExperimentalTorsions (patched)
chiral_sets / tetrahedral_centers: [{ center; atoms[4]; volLo, volHi; fused_small_rings }]
impropers:     [{ atoms[4]; C0, C1, C2, fc }]                       — UFF inversion
bonds, angles                                                       — Stage-C 1-2/1-3 constraints
double_bond_ends, stereo_double_bonds                               — acceptance checks
pin_atoms, core_seeds, sidechain_confs                              — the core-pin recipe
```

## The engine (`bb-embed`)

Per conformer, the faithful ETKDG staged distance-geometry minimization:

- **Init** — metric-matrix eigenvalue embedding (RDKit `useRandomCoords=false`) for core seeds,
  random-box for sidechains. This sets the ensemble spread (the largest fidelity lever we found).
- **Stage A** (4D) — distance-bound violation + chiral volume + 4th-dimension penalty, minimized with
  L-BFGS.
- **Stage B** — 4th-dimension squeeze (chiral weight 0.2, fourth 1.0).
- **Stage C** (3D) — M6 experimental torsions + UFF impropers + distance/angle constraints (flat-bottom),
  L-BFGS.
- **Acceptance checks + retry** — tetrahedral, chiral, planarity, double-bond linear/stereo, final
  chiral bounds; a failing conformer is regenerated, matching RDKit's `embedPoints` retry loop.

The **two-stage core-pin recipe** (matching `build_ligands.py`): `EmbedMolecule` core seeds → freeze the
macrocycle ring via a coordMap (+ re-smoothed bounds) → `EmbedMultipleConfs` sidechains. This clusters
the ensemble to the oracle's diversity.

Numerics delegate to best-in-class crates (argmin L-BFGS, nalgebra eigendecomposition, rand); the ETKDG
force-field terms, staged sequence, checks, and recipe are the reimplementation. Every term is matched
to RDKit source and finite-difference-verified; the triangle-smoothing port is bit-exact vs RDKit; a
force-field-eval cache makes the L-BFGS line search's double-evaluation free.

## Why CPU (the GPU is archived)

The project began targeting the GPU. A full CubeCL→Metal GPU pipeline was built and passes the guard,
but end-to-end it *loses* to the tuned CPU engine: the GPU path is retry-bound (FIRE-minimized
geometries fail the acceptance checks more often, so the retry loop dominates), while the CPU's high
first-attempt acceptance and zero launch overhead win. A residual intra-kernel race (rare,
aggregate-harmless, needs CUDA `racecheck` to pin down) also remains. Several CPU accelerations
(SoA-over-conformers SIMD, spatial active-set, f32 force field) were prototyped and *measured not to pay
off*; the one free win — `lto=fat` + `target-cpu` build flags (~1.3×) — is committed and guard-validated.
The GPU implementation and prototypes live in `.archive/` (gitignored).

## Workspace layout

```
rust/
  bb-core/    MoleculeSpec data contract (serde)                     — no deps beyond serde
  bb-rdkit/   setup: SMILES → MoleculeSpec via cxx FFI (patched RDKit); panic-free library
  bb-embed/   the engine: force field, staged minimize, checks, core-pin recipe; CLIs bb-embed/bb-batch
  bb-py/      thin PyO3 bindings → the `betterbuilder` Python module (maturin)
```

## The Python API (`bb-py`)

A maximally thin PyO3 wrapper — `betterbuilder.embed(smiles, seed) -> Conformers` — with all logic in
Rust (`build_spec` + `embed_recipe`). Invalid SMILES raise `ValueError`. This is the vehicle for the
eventual DOCK-pipeline drop-in (SDF/db2 output; the remaining "output gap").

## Boundaries

- **Perception:** RDKit at the boundary (via the patched-RDKit FFI); native perception is out of scope.
- **Algorithm:** reproduced, not improved — validated by ensemble RMSD, not bit-exact (we use the
  correct math where RDKit has numerical quirks, e.g. dropped factors-of-two now matched).
- **Solvation/charges:** external, out of scope. **Local-validity quality work:** parked (`docs/01`).

## References

- Sampling algorithm: **ETKDGv3** (Wang, Witek, Landrum, Riniker, *JCIM* 2020).
- Kernel math source: RDKit DistGeom / CrystalFF / DistGeomHelpers (the sole oracle; never transcribed
  as a shortcut — the real code is the source of truth).
- Throughput premise: macrocycle build is generation-limited (Foloppe et al., *Bioorg. Med. Chem.* 2019).
