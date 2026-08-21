# db2 format gate — reference corpus

`reference/*.db2` are **golden db2 files from Divya's own pipeline**, used to prove the `bb-db2`
writer is byte-exact:

```
rust/target/release/bb-db2-roundtrip validation/db2/reference/
# db2 round-trip: 9 byte-identical, 0 differ, 0 unreadable (9 files)
```

The gate reads each file, parses it, re-serializes it, and requires the output to be **byte-identical**
to the input. Since these are real oracle output, a byte-identical round-trip certifies `bb-db2`
reproduces the current `mol2db2` writer exactly. This is the *format* layer only; building a
`Db2Entry` from BetterBuilder's own conformers (the `mol2db2` hierarchy) is a separate, later stage.

## Provenance (reproducible)

Generated with the current oracle inside the pipeline container image (read-only reference):
`/nfs/home/amurray2/toolchain/images/build_macrocycle_final.sif`, `$DOCKBASE=/soft/DOCK-macro-latest6`.

1. `build_ligands.py` (Divya's driver) run on the SMILES, producing per-molecule
   `solv/<idx>/output.{mol2,solv}` (RDKit ETKDGv3 embed → obabel → AMSOL → `process_amsol_mol2`).
2. The db2 writer run directly on those, isolating the mol2db2 stage:
   ```
   python $DOCKBASE/ligand/mol2db2_py3_strain/mol2db2.py \
       -m output.mol2 -s output.solv \
       -d $DOCKBASE/ligand/mol2db2/clashfile.txt -o out.db2gz
   ```
   (`mol2db2.py`'s `-o` output is gzip; decompressed here to text. The driver's own
   `mol2db2_quick` returns the same content as plain text.)

`mc0001`–`mc0008` are the first eight macrocycles of `validation/seeds_100.smi`; `benzoic.0` is
`OC(=O)c1ccccc1`. The writer authority is `hierarchy.py` `_allButSetWriter`/`_setWriter` in that image.

Note: `data/zacks_data/*.db2` (3287 files) are real DOCK output from an **older** mol2db2 whose M1
`name`/`protName` are left-justified; the current writer right-justifies them (`%16s`/`%9s`), so those
are not byte-references for the current format.
