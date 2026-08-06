# BetterBuilder — local-validity methodology

Why the validity gate is layered the way it is, grounded in experiments rather than assertion. The
gate flags the deterministically-checkable "local failure modes" of macrocycle conformers — the
class fixable without changing the sampling algorithm — using three layers:

1. **PoseBusters** — general small-molecule geometry (bond lengths/angles vs distance-geometry
   bounds, aromatic & sp2 planarity, internal clash, stereo-consistency, energy). Off the shelf.
2. **CSD torsion library** — the Shoichet STRAIN filter (`Torsion_Strain` + `TL_2.1_VERSION_6.xml`),
   used unmodified.
3. **Native ω / stereo / clash** — peptide-aware backbone checks implemented here.

The obvious question is whether layers 1 and 2 make layer 3 redundant. Two experiments show they do
not: each established tool has a **structural blind spot** exactly where the macrocycle-peptide
pathology lives. Both experiments are reproducible (`experiments/`), on the vendored tools.

## Experiment 1 — PoseBusters does not check amide torsion

Force a cis (ω=0°) and a twisted (ω=90°) amide into a valid peptide conformer, leaving bonds and
angles untouched and inducing no clash, then compare PoseBusters vs the native ω check
(`experiments/pb_vs_omega.py`):

| case          | ω (deg) | native cisNP/twist/clash | PoseBusters |
|---------------|--------:|--------------------------|-------------|
| trans (valid) |   167.6 | 0 / 0 / 0                | all pass    |
| forced CIS    |     0.0 | **1** / 0 / 0            | **all pass**|
| forced TWIST  |    90.0 | 0 / **1** / 0            | **all pass**|

A spurious cis peptide bond, and a 90° twisted amide, both **pass PoseBusters completely.** The
reason is structural: PoseBusters checks 2-body (lengths), 3-body (angles), specific-group planarity,
and clashes. Amide cis/trans/twist is a **4-body torsion** around the C–N bond — it changes no length,
no angle, and needn't clash. PoseBusters deliberately does not encode "peptide bonds should be trans"
(that is peptide domain knowledge, outside a general validator's scope).

## Experiment 2 — the torsion library does not score ring torsions

The knowledge-based torsion library (Guba/Rarey; Shoichet STRAIN) scores torsional strain from CSD
angle histograms. Two facts emerge (`experiments/torsion_library_coverage.py`):

**Part A — acyclic amide (linear peptide):** it flags a twisted amide but barely a cis one.

| case  | ω (deg) | tE (total) | pE (worst torsion) |
|-------|--------:|-----------:|-------------------:|
| trans |   167.6 |       3.91 |               1.24 |
| CIS   |     0.0 |       4.67 |          **2.00**  |
| TWIST |    90.0 |       8.26 |          **5.60**  |

The twist is clearly strained (5.60 vs 1.24 baseline); the cis is barely above trans (2.00 vs 1.24),
because **cis amides are populated in the CSD** — they are valid low-strain minima, just chemically
anomalous for a non-proline residue. (For contrast, an MMFF force field scores this cis at ~6.9
kcal/mol — it would falsely suggest the library flags cis. That disagreement is exactly why a force
field is not a substitute for the knowledge-based library.)

**Part B — the same amide as a ring bond:** the library is blind to it.

| molecule                    | rot. bonds | amides in ring | strain tE | pE |
|-----------------------------|-----------:|---------------:|----------:|----:|
| linear (acyclic amides)     |          8 |            0/3 |      3.91 | 1.24 |
| cyclo(L-Ala)6 (ring amides) |          6 |            6/6 |  **0.00** | **0.00** |

The identical amide reads **3.91 TEU when acyclic and 0.00 when in the ring.** The library scores
only **acyclic rotatable bonds**; a cyclic peptide's entire backbone — every amide ω — is *ring*
bonds it does not score. On real macrocycle conformers it returns ~0 strain across the board (the
gate's own demo: 0/90 conformers flagged), blind to backbone twists the ω check flags ~14% of the
time. This is not a version gap or a bug — it is what a rotatable-bond torsion library does.

## Settled coverage map (macrocyclic peptides)

| failure mode | PoseBusters | Torsion library / STRAIN | native ω (backbone-aware) |
|---|:--:|:--:|:--:|
| bond / angle / planarity / clash | ✅ | — | — |
| exocyclic / side-chain unusual torsion | ✗ | ✅ | ✗ |
| backbone **twisted** amide (in ring) | ✗ | ✗ (blind to ring) | ✅ |
| backbone **cis** amide (in ring) | ✗ | ✗ (cis is CSD-populated) | ✅ |

## Design conclusion

- **The native backbone-ω check is load-bearing, not a redundant hand-roll.** It covers the exact
  gap both PoseBusters and the torsion library structurally leave open for macrocycle rings, and is
  the only source for the cis and twisted *backbone* amides. For the Rust port, `bb-validate`'s
  dihedral measurement over the ring backbone is the essential, must-build component.
- **PoseBusters** stays as the general-geometry backbone; **the torsion library** stays as a bolt-on
  for exocyclic/side-chain torsions (honestly labeled "acyclic torsions only" in the gate output).
- **MolProbity** (Richardson lab) *would* flag cis-nonproline and twisted peptides with
  literature-calibrated frequencies, but is protein-residue-oriented and unreliable on non-canonical
  macrocycles (N-methylation, D-residues, depsipeptide linkers) — a possible cross-check for
  peptide-perceivable cases, not a general macrocycle solution.
- **Calibration:** the cis check is unambiguous (ω near 0 at a non-proline N). The twist threshold
  (currently >30° off planar) is the one knob to calibrate against known-good macrocycle structures
  (e.g. CREMP CREST/xTB ensembles).

## References

- **PoseBusters** — Buttenschoen, Morris, Deane, *Chem. Sci.* 2024, DOI 10.1039/D3SC04185A.
- **Torsion Library (Reloaded)** — Guba, Meyder, Rarey, Hert, *J. Chem. Inf. Model.* 2016, PMID 26679290.
- **STRAIN filter** — Gu, Smith, Bender, Krylov, Coleman, Shoichet, *J. Chem. Inf. Model.* 2021,
  DOI 10.1021/acs.jcim.1c00368 (the vendored `Torsion_Strain` + `TL_2.1_VERSION_6.xml`).
- **MolProbity** — Williams et al., *Protein Sci.* 2018, PMID 29067766.
- **Macrocycle docking is generation-limited** — Foloppe et al., *Bioorg. Med. Chem.* 2019, PMID 31771798.
- **CREMP** (macrocyclic-peptide reference ensembles) — Grambow et al., *Sci. Data* 2024, DOI 10.1038/s41597-024-03698-y.
- **ETKDGv3** (the reproduced sampling algorithm) — Wang, Witek, Landrum, Riniker, *J. Chem. Inf. Model.* 2020, DOI 10.1021/acs.jcim.0c00025.
