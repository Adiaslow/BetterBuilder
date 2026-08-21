#!/usr/bin/env python3
# Corpus-scale db2 WRITER gate (goal criterion 1: output writer byte-identical across the corpus).
# Stage-isolated + per-block, mirroring build_ligands.py: for each molecule, `bb-emit-inputs` produces
# native's per-core-seed-block (mol2, native db2) at file precision (%.4f coords, %8.2f solv, strain
# header); her container `mol2db2.py` writes the golden db2 from the SAME per-block mol2+solv; we
# byte-diff. A block PASSES iff the ONLY difference is the 2 SMILES/longname M-records (which a mol2
# cannot carry — her mol2db2 defaults them to "fake"; native's db2 has the real ones and is correct).
#
# Host (AMSOL via bb-emit-inputs) + container (mol2db2). RDKit-free on the native side.
#   python3 validation/db2/corpus_writer_gate.py [corpus.smi] [N] [workroot]
import gzip
import os
import subprocess
import sys
import pathlib

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "validation/lib"))
import env as oracle_env  # single source: container path (IMG) + AMSOL toolchain
EMIT = REPO / "rust/target/release/bb-emit-inputs"


def is_metadata_only(golden: bytes, native: bytes) -> bool:
    """True iff golden==native except for the 2 SMILES/longname M-records (the mol2-can't-carry-SMILES
    artifact). Any other difference (or a length mismatch) is a real writer divergence → not metadata."""
    gl = golden.decode("latin1").splitlines()
    nl = native.decode("latin1").splitlines()
    if len(gl) != len(nl):
        return False
    diffs = [(a, b) for a, b in zip(gl, nl) if a != b]
    if not diffs:
        return True
    # the only allowed diffs: M-records where golden is the "fake" placeholder
    return all(a[:1] == "M" and a.split()[1:2] == ["fake"] for a, b in diffs)


def run_mol2db2(work: pathlib.Path, name: str, nblocks: int):
    """Run her container mol2db2.py on each block's mol2 → golden .db2gz."""
    script = (
        "cd %s; export PYTHONPATH=$DOCKBASE/ligand/mol2db2_py3_strain; "
        "for b in $(seq 0 %d); do "
        "python3 $DOCKBASE/ligand/mol2db2_py3_strain/mol2db2.py -m %s.b$b.mol2 -s %s.solv "
        "-d $DOCKBASE/ligand/mol2db2/clashfile.txt -o %s.b$b.golden.db2gz 2>/dev/null; done"
        % (work, nblocks - 1, name, name, name)
    )
    subprocess.run(oracle_env.apptainer_exec(["sh", "-c", script], binds=[f"{REPO}:{REPO}", f"{work}:{work}"]),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    corpus = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "validation/seeds_100.smi"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    workroot = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else REPO / "validation/db2/.corpus_gate"
    workroot.mkdir(parents=True, exist_ok=True)

    mols = []
    for line in open(corpus):
        t = line.split()
        if t:
            mols.append((t[0], t[1] if len(t) > 1 else "mol%d" % len(mols)))
        if len(mols) >= n:
            break

    env = oracle_env.amsol_env()
    tot_blocks = tot_ident = tot_meta = tot_diff = 0
    mol_pass = 0
    for smi, name in mols:
        work = workroot / name
        if work.exists():
            for f in work.glob("*"):
                f.unlink() if f.is_file() else None
        work.mkdir(parents=True, exist_ok=True)
        r = subprocess.run([str(EMIT), smi, name, str(work), "210185"], env=env,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True)
        if r.returncode != 0:
            print(f"  {name:16s} EMIT-FAIL ({r.stderr.strip()[:60]})")
            continue
        blocks = sorted(work.glob(f"{name}.b*.native.db2"))
        if not blocks:
            print(f"  {name:16s} no blocks")
            continue
        run_mol2db2(work, name, len(blocks))
        ident = meta = diff = 0
        for nb in blocks:
            b = nb.name.split(".b")[1].split(".")[0]
            g = work / f"{name}.b{b}.golden.db2gz"
            if not g.exists():
                diff += 1
                continue
            try:
                golden = gzip.decompress(g.read_bytes())
            except Exception:
                golden = g.read_bytes()
            native = nb.read_bytes()
            if golden == native:
                ident += 1
            elif is_metadata_only(golden, native):
                meta += 1
            else:
                diff += 1
        nb_tot = len(blocks)
        ok = (diff == 0)
        mol_pass += ok
        tot_blocks += nb_tot; tot_ident += ident; tot_meta += meta; tot_diff += diff
        print(f"  {name:16s} {nb_tot:2d} blocks: {ident:2d} exact, {meta:2d} meta-only, {diff:2d} DIVERGE "
              f"{'PASS' if ok else 'FAIL'}")

    print(f"\n=== corpus writer gate: {mol_pass}/{len(mols)} molecules pass "
          f"(all blocks byte-identical modulo SMILES metadata) ===")
    print(f"  blocks: {tot_ident} fully byte-identical + {tot_meta} identical-except-SMILES "
          f"+ {tot_diff} real divergences  / {tot_blocks} total")
    return 1 if tot_diff else 0


if __name__ == "__main__":
    sys.exit(main())
