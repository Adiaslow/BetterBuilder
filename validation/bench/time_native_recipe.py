#!/usr/bin/env python3
# Speed benchmark (goal dimension 4): time native's conformer-generation path
# (bb-spec-native | bb-embed recipe) per molecule, host-side, RDKit-free. Pairs with Divya's
# per-molecule times (parsed from the gen_reference log / her ensembles) to report the speedup.
#
# THREADING FAIRNESS: native's embed_recipe uses rayon (multi-thread per molecule); her RDKit embed is
# single-threaded (EmbedMolecule/EmbedMultipleConfs, numThreads default 1). Comparing native-multi to
# her-single measures native's PARALLELISM (~20-30x on this 128-core box), NOT algorithmic speed, and
# inflates the ratio ~20x. For the honest per-core / whole-corpus-throughput speedup, run native
# SINGLE-threaded (BB_SINGLE_THREAD=1 → RAYON_NUM_THREADS=1) so both sides use one core. The multi-thread
# time is still recorded (single-molecule latency, a real but separate capability her code could also
# get via numThreads).
#
#   validation/bench/time_native_recipe.py <corpus.smi> [out.json]   [BB_SINGLE_THREAD=1 for the fair number]
import json
import os
import statistics
import subprocess
import sys
import time
import pathlib

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
SPEC = REPO / "rust/target/release/bb-spec-native"
EMBED = REPO / "rust/target/release/bb-embed"


def main():
    corpus = sys.argv[1] if len(sys.argv) > 1 else str(REPO / "validation/seeds_100.smi")
    out = sys.argv[2] if len(sys.argv) > 2 else str(REPO / "validation/bench/native_times.json")
    single = os.environ.get("BB_SINGLE_THREAD") == "1"
    env = dict(os.environ, RAYON_NUM_THREADS="1") if single else None
    print(f"native timing: {'SINGLE-thread (fair per-core vs her single-thread)' if single else 'MULTI-thread (rayon; latency, NOT a fair speedup vs her single-thread)'}")
    rows = []
    for line in open(corpus):
        tok = line.split()
        if not tok:
            continue
        smi = tok[0]
        name = tok[1] if len(tok) > 1 else "mol%d" % len(rows)
        t0 = time.time()
        spec = subprocess.run([str(SPEC), smi], stdout=subprocess.PIPE).stdout
        conf = subprocess.run([str(EMBED), "-", "recipe", "210185"], input=spec,
                              stdout=subprocess.PIPE, env=env).stdout
        dt = time.time() - t0
        try:
            d = json.loads(conf)
            n_conf = len(d["conformers"])
            n_atoms = d["n_atoms"]
        except Exception:
            n_conf, n_atoms = 0, 0
        rows.append({"name": name, "seconds": round(dt, 3), "n_conformers": n_conf, "n_atoms": n_atoms})
        print("  %-14s %6.3fs  %3d confs (n=%d)" % (name, dt, n_conf, n_atoms))
        sys.stdout.flush()
    json.dump(rows, open(out, "w"), indent=1)
    ts = sorted(r["seconds"] for r in rows)
    med = statistics.median(ts) if ts else 0.0
    print("native recipe: %d molecules, median %.3fs, max %.3fs (%s)  -> %s"
          % (len(rows), med, max(ts) if ts else 0, next((r["name"] for r in rows if r["seconds"] == max(ts)), ""), out))


if __name__ == "__main__":
    sys.exit(main())
