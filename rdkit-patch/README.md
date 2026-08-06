# Patched RDKit — the `BB_RDKIT_ROOT` build dependency

BetterBuilder's setup layer (`rust/bb-rdkit`) links a **patched RDKit build**, not stock RDKit. This
directory is the patch that defines it — so a clean checkout can reproduce the dependency that
`BB_RDKIT_ROOT` points at.

## What the patch does

It raises the **macrocycle amide-bond torsion force constants** in RDKit's CrystalFF experimental-
torsion table so amide bonds stay planar during ETKDG embedding (stock RDKit under-penalizes this for
macrocycles). Concretely, in `torsionPreferences_macrocycles.in`:

- ring amide patterns `[C][C;r{9-}](=O)@;-[NX3H*;r][CX4H*]` : force constant **8.0 → 80.0**
- the `[O]=[CX3;r{9-}](a)@;-[NX3H0;r][!#1]` amide-planarity pattern : **13.9 → 130.9** (and added)

This is why an extracted `MoleculeSpec` carries amide torsion FC 80.0 (asserted by
`bb-rdkit`'s `spec_carries_torsions_and_chirality` test) rather than stock 8.0.

## Contents

| File | What it is |
|---|---|
| `divya_amide_torsions.diff` | The change as a git diff against RDKit `Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionPreferences_macrocycles.in`. |
| `torsionPreferences_macrocycles.in` | The full **patched** table (drop-in replacement). |
| `torsionPreferences_v2.in` | The companion chain-torsion table (unchanged upstream copy, pinned for the build). |

## Building the dependency

Against an RDKit source tree (BetterBuilder targets **RDKit 2026.09.1pre**):

```sh
RDKIT_SRC=/path/to/rdkit           # a git checkout of RDKit at the target version
CF=$RDKIT_SRC/Code/GraphMol/ForceFieldHelpers/CrystalFF

# apply the patch (either the diff, or drop in the patched tables):
cp torsionPreferences_macrocycles.in torsionPreferences_v2.in "$CF"/
# or:  git -C "$RDKIT_SRC" apply /path/to/divya_amide_torsions.diff

# build RDKit from source (see RDKit's INSTALL docs), then point BetterBuilder at it:
export BB_RDKIT_ROOT="$RDKIT_SRC"  # the dir containing Code/ and lib/
```

Then build BetterBuilder normally (`cd rust && cargo build --release`). See the root `README.md`.
