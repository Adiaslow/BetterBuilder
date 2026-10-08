# Parity harness

Checks that BetterBuilder reproduces the oracle's behaviour, stage by stage, upstream to downstream.
BetterBuilder replaces the oracle's Python with Rust; the machinery differs deliberately. What must
match is the output behaviour.

## The oracle

Divya's unmodified pipeline, run in her **production** container, `build_macrocycle_final.sif` —
the image her Wynton job scripts invoke (`apptainer exec ... --bind ${INDIR}:/data ...
build_macrocycle_final.sif`), confirmed with her directly.

Do not pick an image by timestamp. `build_mc_july_15.sif` is newer and its `build_ligands.py` matches
the `DOCK-macro-latest6` tree on the shared filesystem byte-for-byte — but that tree is her working
copy, and the July image is a variant that uses corina rather than RDKit for the solvation starting
conformer. The two images share a byte-identical conformer-generation block (the 138 lines from
`count_exo_rotatable` to `###END RDKIT`) and differ upstream of it.

Reference values come only from her code executing. Three techniques, none of which reimplement it:

| Technique | Used for | How |
|---|---|---|
| **Direct library calls** | RDKit perception, `mol2`/`solv`/`hierarchy` state | import her modules, call them, read the results |
| **AST extraction** | `count_exo_rotatable` | locate the `FunctionDef` in her file, compile and execute that node — her exact bytes |
| **Frame tracing** | `cmap` (pinned atoms), conformer ensembles | run her file unmodified under `sys.settrace`, read locals from the live frame |

Tracing is scoped to her file (`return None` for other frames); tracing into RDKit slows a run by
orders of magnitude.

Nothing BetterBuilder produces is ever a reference, and neither is any reimplementation of hers —
including one written for the purpose. Values computed by a modified copy of the oracle are not
oracle values.

## Contents

| File | Role |
|---|---|
| `dump_perception_batch.py` | RDKit's perception of many molecules (stage 1) |
| `dump_oracle.py` | `mol2`/`solv`/`hierarchy` internal state and db2 bytes (stages 4-6) |
| `dump_valence.py` | RDKit's per-element valence lists, from its own accessors (implicit hydrogens) |
| `dump_kekulize.py` | RDKit's kekulized bond orders and perceived aromaticity |
| `dump_addhs.py` | RDKit's molecule after `AddHs`: atom and bond order as sequences |
| `dump_hybrid.py` | RDKit's conjugation and hybridization |
| `dump_hcount.py` | RDKit's per-atom valences and hydrogen counts over a whole corpus |
| `dump_rings.py` | RDKit's symmetrized SSSR over a whole corpus (ring perception) |
| `dump_exo.py` | her `count_exo_rotatable` extracted by AST and run over a whole corpus (stage 2) |
| `trace_pins.py` | the coordMap her seed loop builds, and the recipe parameters in the same frame (stage 2) |
| `trace_confs.py` | her conformer ensembles (stage 3) |
| `gate_spec.py` | native `MoleculeSpec` (`bb-spec-native`) vs RDKit-bridge spec (`bb-spec`), field by field |
| `dump.sh` | runs `dump_oracle.py` over a completed pipeline run |
| `fixtures/` | captured oracle values |

## The corpus

A gate reports only on the chemistry it was run against, so the corpus is an explicit input. Each
selector is a glob, so coverage extends by adding a provenanced corpus and its fixture rather than
by editing `check.py`: `BB_CORPUS` (`seeds_*.smi`), `BB_EXO` (`exo_*.json`), `BB_STAGE2`
(`stage2_*.json`), `BB_PERCEPTION`, `BB_ENSEMBLES`.

`validation/seeds_5000.smi` is a systematic 1-in-1800 sample of the production workload,
`/wynton/group/bks/work/kholland/resampled_synthMC_fp_lib/9M_macrocycle_seeds.smi` (kholland, bks,
read-only). `seeds_5000.provenance` records source, size, sampling method, date and md5. Everything
else is drawn from it: `seeds_100.smi`, the traced 25, and `seeds_branch20x10.smi` — the last
selected to reach the `exo<=2` recipe branch, which 295 of the 5000 take.

The property that matters for this pipeline is Divya's amide reweighting, since it is the only
reason her RDKit is not stock. Its four SMARTS fire on **83%** of the production workload, at 1.75
patched torsions per molecule. Any corpus is characterised by that number before it is used.

## Running

Her production `build_ligands.py` reads protomers from **stdin** as `<smiles> <name>` and assigns
protomer ids itself; it also loads `/data/input.smi`, a neutral-SMILES dictionary requiring unique
names. (The `sys.argv[1]` three-column path in that file is commented out — it is the older
contract, and the one the `DOCK-macro-latest6` tree copy still uses.)

```sh
apptainer exec --pwd /work/out -B <run>:/work -B <run>/data:/data -B validation/parity:/parity \
    <image> python3 /parity/trace_confs.py /work/in/protomers2col.smi /parity/fixtures/ensembles_seeds
```

It skips solvation when `solv/<mol>/` already exists, and ends by moving `output.tar.gz` to `/data`
— bind something there or expect the last step to fail after the work is done.

Stage 4 needs a completed run plus AMSOL. AMSOL links `libg2c.so.0`, which is not on the host, so
extract it and its libraries from the image once and point `bb-solv` at them:

```sh
BB_ORACLE_RUN=<run> BB_AMSOL_EXE=<dir>/amsol7.1 BB_AMSOL_LD_LIBRARY_PATH=<dir>/lib \
    validation/parity/check.py stage4
```

The extracted binary was checked against hers by running her own `temp.in-wat` through it: output
byte-identical apart from the wall-clock line.

`trace_confs.py` locates its capture point by source text (`tmp_out_sdf = tempfile.NamedTemporaryFile`)
rather than a line number, and matches her file by resolved path or basename: `runpy` reports a
`co_filename` that does not string-equal the literal path, and an exact comparison rejects every
frame and captures nothing without erroring.

## Status

Measured on the production seed sample: `seeds_5000.smi` for the recipe check, `seeds_100.smi` (83%
patch coverage) for perception, and traced subsets of 25 and 12 (88% patch coverage) for pins and
ensembles.

Rows 0–4 were measured before the 2026-10 changes. Stage 3 needs re-measurement: the recipe now
replaces a core that cannot hold its sidechains and gives each embed RDKit's attempt budget (10 × atom
count). The one perception change of that round (no double-bond stereo in rings under 8 atoms, as
RDKit) disagrees with the previous behaviour on none of the 6,000 corpus molecules, which have no such
ring. Rows 5 and 6 were run on 2026-10-07.

| Stage | Scope | Result |
|---|---|---|
| 0. Ring perception (Rust) | 5000 molecules | **exact** — every ring identical to RDKit's symmetrized SSSR |
| 0. Valence & hydrogens (Rust) | 5000 molecules, 205683 atoms | **exact** — explicit valence, implicit valence and total H, per atom |
| 0. AddHs (Rust) | 5000 molecules, 403420 atoms | **exact** — atom order and bond order, both as sequences |
| 0. Aromaticity (assumption) | 5000 molecules, 41547 aromatic atoms | written aromaticity equals RDKit's perceived, 0 mismatches |
| 0. Recipe & pinned core (Rust) | 5000 molecules / 36 traced | **exact** — `count_exo_rotatable`, `core_seeds`/`sidechain_confs`, largest ring |
| 0. Angles (Rust) | 600 molecules, 94668 angles | **exact** — sequence identical, including the linear flag |
| 0. Conjugation & hybridization (Rust) | 5000 molecules, 403420 atoms | **exact** — per atom and per bond |
| 1. Perception | 100 molecules | exact — bounds within 2e-6 A (f32 storage); every RDKit library torsion present with identical atoms and force constants |
| 2. Recipe parameters | 5000 molecules, both branches | exact — 0 mismatches against her extracted `count_exo_rotatable` |
| 2. Pinned atoms | 36 molecules, both branches | exact — pin sets identical, and constant across each molecule's core seeds |
| 3. Conformer ensembles | 23 of 25 molecules | NN 1.03x the oracle's own seed-to-seed spread; within-seed diversity 1.04x, between-seed 1.08x |
| 4. Solvation | 24 molecules | **byte-identical** `.solv` given her geometry |
| 5. Strain | 9 molecules | total and max strain within 1e-3 of her `Torsion_Strain`, given her geometry (`bb-strain` gate) |
| 6. db2 | 13 references | **byte-identical** db2 and identical hierarchy internals given her molecule and solvation; mol2 byte-identical on the 9 single-conformation ones (`bb-output` gate) |

Ring perception is the first piece of the Rust front-end (`bb-perceive`, parsing via `yowl`), and it
reproduces RDKit's symmetrized SSSR exactly on the corpus. Getting there required the traversal, not
a closed-form rule: degree-2 chains contribute one search root each, rings are kept the first time
their atom set is seen, roots that discover the same ring are re-searched with each other's bonds
removed, and the survivors are pruned by bond coverage.

**Bond order is part of that traversal, not a detail.** Ring search walks an atom's neighbours in
bond order, so among equal-size rings it decides which is found. RDKit's order comes from its parser:
`SmilesParse.cpp:256` calls `CloseMolRings` only after the whole string is parsed, so every chain
bond precedes every ring-closure bond, and `CloseMolRings` walks atom bookmarks keyed by
ring-closure digit, so closures are ordered by digit rather than by position in the string.
Reproducing that ordering took the corpus from 99.96% to exact. A bond list that matches as a *set*
is not enough; it has to match as a sequence.

That sensitivity is worth keeping in mind downstream: swapping RDKit's plain SSSR for its
symmetrized one moves the distance bounds by up to 1.9 A, and a single differing ring moves them by
~0.06 A — against stage 1's 2e-6 A tolerance. The bounds matrix is therefore a far sharper detector
of a ring or ordering error than ring-set equality is.

A minimum-cycle-basis implementation is kept in `rings.rs` as an independent cross-check; on its own
it diverges on `r{9-}` for 7 of 5000, because a basis omits the symmetry-equivalent rings that
predicate needs.

Stage 1 compares `bb-spec` against RDKit's Python API, and `bb-spec` reaches RDKit through the cxx
bridge — so it establishes that the bridge forwards RDKit's values, including the patched amide
torsions, on chemistry where they fire. It is not a fidelity gate for the Rust reimplementation of
perception; that is the stage-0 rows and `rust/gates.sh`.

Stages 2 and 3 read values out of live frames of her `build_ligands.py` executing. Stage 1 cannot:
the bounds matrix is built inside `EmbedMolecule` in C++ and never reaches a Python frame. It is
therefore her patched RDKit, in her production container, under parameters this harness sets — and
those parameters are equivalent to hers rather than taken from her. `ETKDGv3()` already defaults
`useMacrocycleTorsions`, `useMacrocycle14config` and `useExpTorsionAnglePrefs` to `True`, so her
assignments at `build_ligands.py:474-475` are redundant, and her sidechain embed at `:493` omitting
`useMacrocycleTorsions` does not change its configuration. Re-check this if either side changes.

Stage 4 feeds `bb-solv` the geometry she actually gave AMSOL (`3d/<name>.mol2`, the file her run
moves out of `solv/`) rather than a conformer of ours. Her solvation embed at `build_ligands.py:191`
is unseeded, so an end-to-end comparison would be distributional; fixing the geometry makes the
chain deterministic and attributes any difference to us. What it therefore covers is Z-matrix
construction, AMSOL input generation, invocation and output parsing — not the choice of conformer
going in, which is stage 4b and is not yet measured. It also runs without OpenEye, which her
`make_amsol71_input_docker.py3.py` needs for `OENetCharge`; `bb-solv` takes formal charge from
perception instead.

Of the 25 molecules traced, one produced 189 of 200 conformers — the conformer-to-seed mapping is
not recoverable from a partial ensemble, so the seed-structured metrics are skipped and the count
reported. One further molecule reached solvation but never conformer generation.

## Choosing a corpus

A gate is bounded by its input, and an unprovenanced input bounds it invisibly. Stages 1-3 were
first measured on a committed molecule set with no recorded origin — a sample of a generated
library, on which the amide patch fired at 16%, against a production workload where it fires at
83%. The oracle really did run her code; the corpus simply never asked it about the chemistry that
distinguishes her engine. The lesson is the same one already applied to reference *values*: record
where the data came from, and characterise it against the property under test before trusting a
result measured on it.

Selecting which molecules to test is part of the corpus, not separate from it. Picking candidates
with our own code samples only where our code already agrees with itself: a disagreement it produces
will surface, but a molecule it wrongly fails to nominate never gets looked at. Where her code can
do the selecting — `count_exo_rotatable` extracted and executed — it should.

The program is a separate axis, and here it happened not to matter: the 138-line
conformer-generation block is byte-identical between the two images, and it starts from
`Chem.MolFromSmiles(mol.smiles)`, so the ensembles do not depend on the corina-vs-RDKit difference
upstream. What does differ is the input contract and every line number, which is why capture points
are located by source text and fixtures are re-derived when the image changes.

## Choosing a metric

Four errors were made and corrected while establishing the above. They share a shape:
comparing two quantities before establishing they are the same kind of thing.

- **The baseline must be the oracle's.** Comparing our disagreement with the oracle against *our own*
  seed-to-seed spread proves nothing: an implementation that sampled too narrowly would show a small
  self-spread and score well. The oracle's own spread is measurable from its output — its 200
  conformers are `core_seeds` blocks of `sidechain_confs`, so first half against second half gives
  its between-seed variation. Match sample sizes; nearest-neighbour distance shrinks as the target
  set grows.
- **Nearest-neighbour distance is one-directional.** It asks whether each of ours is near something
  of theirs, not whether we reproduce their spread. An ensemble collapsed into one basin scores well.
  Check the reverse direction, and compare diversity directly — within-seed (sidechain sampling
  around a pinned core) and between-seed (core conformations). Matching the oracle's spread is the
  requirement; being more homogeneous is a defect, not an achievement.
- **Identify the oracle by what production runs, not by timestamp.** The newest image, and the one
  whose source matches the shared filesystem, was a variant rather than the deployed pipeline. The
  authority is the job script's `apptainer exec` line — or asking.
- **Establish what each side is reporting.** RDKit's Python `GetExperimentalTorsions` returns only
  SMARTS-matched library torsions, while `CrystalFFDetails` — what the bridge reads — also carries
  the basic-knowledge planarity terms. Equal counts is the wrong relation; containment is the right
  one. Filtering by a guessed signature (`V=[0,100,0,0,0,0]`) is wrong too, because some library
  torsions share it.

## Open

- Stage 3 covers 23 molecules, stage 1 covers 100. Neither is the full workload; the stage-2
  recipe check covers 5000.
- Stage 4b — solvation from a conformer of ours rather than hers — is not measured. Because her
  starting embed is unseeded, the baseline must be her own run-to-run spread, which means running
  her solvation twice over the same molecules first.
- The Rust front-end now assembles the whole `MoleculeSpec` with no RDKit: SMILES parsing, ring
  perception, valence/hydrogen counts, `AddHs`, the recipe parameters, the pinned core, bond angles,
  conjugation and hybridization, torsion assignment, chiral sets, stereo double bonds, impropers, and
  the full `setTopolBounds` matrix. `gate_spec.py` runs the native assembler (`bb-spec-native`) and
  the RDKit bridge (`bb-spec`) on each SMILES and compares every field of the emitted spec; the
  bridge is the oracle. It matches field-for-field on all 5000 seeds (0 differ), with `stereo_double_bonds`
  compared up to reference-atom equivalence (see below). `bb-spec-native` links no RDKit or boost —
  `ldd` shows neither — so the spec path is RDKit-free end to end. The bridge remains only as this
  gate's oracle, not in the pipeline.
- `gate_spec.py` at full corpus is what caught two fused-heteroaromatic molecules (seeds 2111, 2894)
  where a bare SMARTS `r<n>` was read as ring membership rather than smallest-ring size — the
  600-molecule torsion fixture did not reach them. RDKit's `r<n>` matches iff an atom's *smallest*
  ring is size n (a fusion atom shared by a 5- and 6-ring fails `r6`); fixing that in `smarts_match.rs`
  brought the corpus to 0 differences. Component fixtures gate a slice; the full spec gate gates the
  corpus.
- Hybridization is assigned to heavy atoms only. RDKit sanitizes before `AddHs`, so appended
  hydrogens keep the default and are unspecified — 197737 of them in this corpus, which is exactly
  the hydrogen count.
- The angles fixture is captured through the bridge, since `collectBondsAndAngles` is not exposed to
  RDKit's Python API. The bridge calls that function and copies its result, so the values are
  RDKit's; it is the one fixture not reachable by a direct library call or by tracing.
- `Bond::IsInRing` means the bond lies on a cycle, which is *not* the same as surviving degree-1
  trimming: a bridge joining two ring systems keeps both atoms at degree 2 and survives, yet lies on
  no cycle. Ring-bond membership is therefore cut-edge detection.
- Kekulization is not implemented and is not needed: 42223 bonds are still `AROMATIC` after
  sanitization, so the molecule reaching `EmbedMolecule` carries aromatic bonds and the Kekule
  assignment never leaves sanitization. Revisit only if something downstream reads integer bond
  orders for aromatic bonds.
- Aromaticity perception is not implemented either; the front-end reads what the SMILES writes. That
  holds for this workload and is gated, not assumed — but it is a property of the corpus, and
  Kekule-form input would break it completely. The gate is the thing that would catch that.
- Hydrogen counts are computed from aromaticity **as written in the SMILES**, which is what
  `updatePropertyCache` sees before kekulization. Re-check them once aromaticity is perceived
  rather than read.
- Stages 5-6 are not measured; strain and db2 are not built.
- `trace_pins.py` and `trace_confs.py` line-trace her whole file to read one line, which costs
  roughly an hour per 25 molecules. Installing a local trace only on the frame that matters would
  make re-capture cheap enough for a larger slice.

## Oracle behaviour observed

Findings about her pipeline, from these runs. Recorded here because they bound what parity means.

- **A molecule is dropped when its hexadecane AMSOL leg exceeds 60 s.** `calc_solvation.py3.csh:10`
  sets `amsoltlimit = 1m` and line 118 runs the leg under `timeout`. On a kill the output is
  truncated at a stdio buffer boundary with no failure marker, and no `output.solv` is written, so
  the molecule vanishes silently. Observed on 1 of 25 production molecules (`CSLB00008uI2Z7`), whose
  water leg completed in 3.38 s.
- **Conformer generation can return fewer than the recipe asks for**: 189 of 200 on one molecule.
- The run opens with `rm: missing operand` — an unguarded shell-out with an empty variable.
