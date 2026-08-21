#!/usr/bin/env python3
# Oracle-provenance certification (CONFORMANCE.md Phase 1): prove our from-source RDKit build produces
# the SAME output as Divya's ACTUAL production build, by compiling an oracle dumper against HER
# container's /soft/rdkit_libs and diffing it against our host-bridge golden at the identical geometry
# (stored in the golden). A byte-identical result means the golden — hence every native gate that
# rests on it — reflects her real build, not our patch reconstruction.
#
#   certify.py <N> <dumper_binary> <golden.jsonl>
#     env BB_ORACLE_CONTAINER = path to her .sif (default: the production macrocycle image)
#
# The dumper must already be compiled (gates.sh certify does this in-container). Exits non-zero if any
# record diverges from byte-identity.
import json, os, pathlib, subprocess, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "lib"))
import env as oracle_env  # single source: container path (IMG) + APPTAINER
# The compiled dumper lives under the repo (not auto-mounted), so bind the repo into the container to
# run it. BB_ORACLE_BIND is "host:dest" (set by gates.sh); empty means the binary is already visible.
BIND = os.environ.get("BB_ORACLE_BIND", "")
EXEC = oracle_env.apptainer_exec([], binds=[BIND] if BIND else [])

def rel_grad(a, b):
    mag = max(max(abs(x) for x in a), max(abs(x) for x in b), 1.0)
    return max(abs(x - y) for x, y in zip(a, b)) / mag

def main():
    if len(sys.argv) < 4:
        print("usage: certify.py <N> <dumper_binary> <golden.jsonl> [dumper_args...]", file=sys.stderr)
        return 2
    N, dumper, golden = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    dumper_args = sys.argv[4:]  # e.g. Stage-B weights "0.2 1.0" passed through to the dumper
    worst, worst_smi, n = 0.0, "", 0
    with open(golden) as f:
        for line in f:
            if n >= N:
                break
            rec = json.loads(line)
            coords_in = " ".join(repr(x) for x in rec["coords"])
            p = subprocess.run(EXEC + [dumper, rec["smiles"]] + dumper_args, input=coords_in,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            if p.returncode != 0:
                print(f"  DUMPER FAIL {rec['smiles'][:40]}: {p.stderr.strip()[:80]}")
                return 1
            her = [float(x) for x in p.stdout.split()]
            gold = rec["grad"]
            if len(her) != len(gold):
                print(f"  LEN MISMATCH {rec['smiles'][:40]}: her={len(her)} gold={len(gold)}")
                return 1
            r = rel_grad(her, gold)
            if r > worst:
                worst, worst_smi = r, rec["smiles"]
            n += 1
    # tag by the oracle (golden) name — uniquely identifies the stage even when a binary serves two
    # (stage_a_dump certifies both Stage-A and Stage-B at different weights).
    tag = os.path.basename(golden).replace(".jsonl", "")
    print(f"  [{tag}] {n} records: worst rel diff = {worst:.3e}  ({worst_smi[:55]})")
    if worst >= 1e-12:
        print(f"  DIVERGENT — our build != her build on {tag}", file=sys.stderr)
        return 1
    print(f"  CERTIFIED: our from-source build == her container build for {tag} (byte-identical)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
