# Patched RDKit — the `BB_RDKIT_ROOT` build dependency

BetterBuilder's setup layer (`rust/bb-rdkit`) links a patched RDKit build. This directory holds the
pin and the patches that define it, so a clean checkout reproduces what `BB_RDKIT_ROOT` points at.

## The pin

Upstream RDKit at commit **`1dfc9b7a1`** (version 2026.09.1pre), from
<https://github.com/rdkit/rdkit>.

## Divergence from upstream

Four files differ, in two groups.

**Torsion tables — these change embedding behavior.** In
`Code/GraphMol/ForceFieldHelpers/CrystalFF/`:

- `torsionPreferences_macrocycles.in`
  - ring amide patterns `[C][C;r{9-}](=O)@;-[NX3H*;r][CX4H*]`: force constant 8.0 → 80.0
  - amide-planarity pattern `[O]=[CX3;r{9-}](a)@;-[NX3H0;r][!#1]`: force constant 13.9 → 130.9
- `torsionPreferences_v2.in`
  - the amide-planarity pattern above is added at force constant 130.9, ahead of the chain-torsion
    entries

An extracted `MoleculeSpec` therefore carries amide torsion force constant 80.0 rather than 8.0,
which `bb-rdkit`'s `spec_carries_torsions_and_chirality` test asserts.

**Compiler compatibility — no behavior change.** These are C++17 spellings of constructs that are
legal in C++20; GCC 13 compiles the tree without them. They are applied anyway so the source matches
the reference build exactly, leaving no judgement call about what is equivalent.

- `Code/GraphMol/QueryOps.cpp`: a structured binding captured by a lambda is replaced with two
  `bool` locals.
- `Code/RDGeneral/JSONHelpers.h`: `static_cast<T::_integral>` gains the `typename` keyword.

## Contents

| File | What it is |
|---|---|
| `divya_amide_torsions.diff` | The torsion-table changes as a git diff; covers both tables. |
| `cxx_compat.diff` | The `QueryOps.cpp` and `JSONHelpers.h` changes as a git diff. |
| `torsionPreferences_macrocycles.in` | The patched macrocycle table, as a drop-in replacement. |
| `torsionPreferences_v2.in` | The patched chain-torsion table, as a drop-in replacement. |

## Building the dependency

`../toolchain/build-rdkit.sh` clones the pin, applies these patches, and builds it;
`../toolchain/build-all.sh` builds the compiler and libraries it needs first. See
[`../toolchain/README.md`](../toolchain/README.md) for prerequisites and options.

The torsion tables can be dropped in instead of applying `divya_amide_torsions.diff`:

```sh
cp torsionPreferences_*.in "$RDKIT_SRC"/Code/GraphMol/ForceFieldHelpers/CrystalFF/
```

`bb-rdkit`'s `spec_carries_torsions_and_chirality` test asserts amide force constant 80.0, so it
fails against an unpatched RDKit.

## `atomic_data.tsv`

RDKit's `periodicTableAtomData` (`Code/GraphMol/atomic_data.cpp`) verbatim: the seven raw-string
chunks concatenated in source order, which is already atomic-number order. Rows are `atomic number,
symbol, row, covalent radius, rB0, van der Waals radius, mass, outer electrons, common isotope,
common isotope mass, valence...`, and a valence of `-1` means unconstrained.

Fields are separated by runs of spaces and tabs with empty tokens dropped, matching
`boost::char_separator<char>(" \t")` in `atomicData::atomicData`. Some rows leave the row-number
column blank, so splitting on single tabs shifts every later column — hydrogen in particular.

`bb-perceive` reads this file directly. `validation/parity/fixtures/valence_rdkit.json`, captured
from RDKit's own accessors, is the check that the parse agrees.
