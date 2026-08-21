"""Summarise an oracle pipeline run: what it produced, what it dropped, and where time went.

    validation/parity/oracle_report.py <run-dir>

<run-dir> is the working directory of a completed `build_ligands.py` run (containing solv/, and
either out/data/output.tar.gz or output.tar.gz), with oracle.log alongside it.

The oracle reports its solvation count but never its final count, and molecules can also be lost
after solvation, so the drop rate is only visible by counting the tarball.
"""
import pathlib
import re
import sys
import tarfile


def find(run, *names):
    for n in names:
        p = run / n
        if p.exists():
            return p
    return None


def main(rundir):
    run = pathlib.Path(rundir)
    log = find(run, "oracle.log", "../oracle.log")
    tar = find(run, "out/data/output.tar.gz", "out/output.tar.gz",
               "data/output.tar.gz", "output.tar.gz")

    delivered, truncated = set(), False
    if tar:
        # The oracle moves output.tar.gz while its writer is still open (build_ligands.py:715,
        # inside the `with tarfile.open(...)`), so when /data is a different filesystem the archive
        # is truncated. Read what is there and say so.
        names = []
        try:
            with tarfile.open(tar, "r|gz") as t:
                for m in t:
                    names.append(m.name)
        except (tarfile.ReadError, EOFError, OSError):
            truncated = True
        for n in names:
            m = re.search(r"(mc\d+)", n)
            if m:
                delivered.add(m.group(1))
        note = "  [TRUNCATED — see build_ligands.py:715]" if truncated else ""
        print(f"  tarball        : {tar.name}, {len(names)} entries, "
              f"{len(delivered)} molecules{note}")

    solvated, failed_amsol = set(), {}
    solvdir = find(run, "solv", "out/solv")
    if solvdir:
        for d in sorted(solvdir.iterdir()):
            m = re.search(r"(mc\d+)", d.name)
            if not m:
                continue
            mid = m.group(1)
            if (d / "output.solv").exists():
                solvated.add(mid)
            else:
                legs = [p for p in d.glob("temp.o-*")
                        if "was not completed successfully" in p.read_text(errors="ignore")]
                failed_amsol[mid] = len(legs)
        print(f"  solvated       : {len(solvated)} of {len(solvated) + len(failed_amsol)}")

    lost_solv = sorted(failed_amsol)
    # A truncated archive loses its tail, so absence from it does not mean the molecule failed.
    # Only conclude post-solvation loss when the archive is intact.
    lost_after = sorted(solvated - delivered) if (delivered and not truncated) else []
    if lost_solv or lost_after:
        print("  dropped:")
        for mid in lost_solv:
            why = "AMSOL size limit" if failed_amsol[mid] else "solvation post-processing"
            print(f"    {mid}: {why}")
        for mid in lost_after:
            print(f"    {mid}: lost after solvation (strain or db2)")
    total = len(solvated) + len(failed_amsol)
    if total:
        lost = len(lost_solv) + len(lost_after)
        print(f"  drop rate      : {lost}/{total} = {100.0 * lost / total:.0f}%")
    if truncated and delivered:
        hidden = sorted(solvated - delivered)
        if hidden:
            print(f"  hidden by truncation: {len(hidden)} molecule(s) built but unreadable "
                  f"({', '.join(hidden)})")

    if log and log.exists():
        text = log.read_text(errors="ignore")
        print("  reported times :")
        for label in ("solvation", "strain", "db2", "rdkit_conf_gen"):
            m = re.search(rf"^{label}:\s+([0-9.]+)", text, re.M)
            if m:
                secs = float(m.group(1))
                per = f"  ({secs / len(delivered):.1f} s/molecule)" if delivered else ""
                print(f"    {label:15s} {secs:9.1f} s{per}")
        print("    note: the oracle's reported `solvation` total is its last molecule's elapsed "
              "time, not the sum (build_ligands.py overwrites the accumulator).")


if __name__ == "__main__":
    main(sys.argv[1])
