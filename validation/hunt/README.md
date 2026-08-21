# Divergence-hunting tier (goal criterion 5)

The corpus gates prove native == her build on **sampled** inputs. This tier makes "no diffs" mean
**"searched hard,"** not "the sample passed" — the CSmith method the conformance contract commits to
(`validation/parity/CONFORMANCE.md` §5, Phase 2). It **generates** adversarial inputs off — but near —
her real workload and diffs native against the RDKit bridge on every one, with **measured coverage**.

## Pieces

- **`mutate_corpus.py`** — structure-aware generator. Seeds from her real macrocycle workload
  (`seeds_5000.smi`) and applies one small, valence-valid, **macrocycle-preserving** mutation per draw
  (heteroatom swap C/N/O/S · substituent add C/N/O/F/Cl · terminal delete · non-ring bond-order flip),
  sanitises, and keeps only valid molecules that still contain a ring ≥ 9. Fixed-seed PRNG →
  deterministic. Prints **measured coverage**: unique count, yield, mutation-op mix, and the ring-size
  / heteroatom-count / atom-count ranges actually reached (so the search breadth is a number, not a
  claim).

- **`differential.py`** — the failable diff. For every generated molecule runs the **bridge**
  (`bb-spec`, RDKit oracle) and **native** (`bb-spec-native`) and compares the `MoleculeSpec`:
  - **HARD** (must match exactly; a diff is a perception bug / spec gap): element list + AddHs order
    (`atomic_numbers`), `bonds`, recipe (`core_seeds`/`sidechain_confs`), `pin_atoms`, and the
    torsion / improper / angle **topology** (the atom-tuple sets — where her amide patch lives).
  - **NUMERIC** (`bounds`, the shared f32 matrix): compared to the vetted tolerance.
  - **STEREO** (`chiral_sets`, `tetrahedral_centers`, `stereo_double_bonds`, `double_bond_ends`):
    reported **separately** because some differ only in representation
    (see the stereo-representation-equivalence note) — a raw mismatch there is a finding to inspect,
    not automatically a defect.
  - **PARSE** divergence: one of {bridge, native} accepts a molecule the other rejects.

  Exits nonzero on any HARD or PARSE divergence, and reports how many molecules were actually diffed.

## Run

```
apptainer exec --bind $REPO:$REPO <build_macrocycle_final.sif> \
  python3 validation/hunt/mutate_corpus.py validation/seeds_5000.smi validation/hunt/generated.smi 2 20210185
apptainer exec --bind $REPO:$REPO <build_macrocycle_final.sif> \
  python3 validation/hunt/differential.py validation/hunt/generated.smi
```

`generated.smi` is a regenerated artifact (deterministic from the seed) — gitignored, not committed,
like the ensemble fixtures. The embed-layer quantities (FF gradients, the six checks, coord_map
bounds, reject energy) are hunted by pointing the existing vetted gate at the same corpus:
`BB_PARITY_CORPUS=validation/hunt/generated.smi cargo test -p bb-rdkit --test corpus_parity`.
