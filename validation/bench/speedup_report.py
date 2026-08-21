#!/usr/bin/env python3
# Goal dimensions 3 (improvement/staller census) + 4 (speedup). Combines:
#   - native_times.json            (native recipe per-molecule wall-clock, host, from time_native_recipe.py)
#   - the gen_reference log         (her per-molecule times + STALLED lines)
#   - ensembles/_stalls.json        (the molecules her pipeline stalls/fails on)
# and reports the per-molecule speedup distribution on clean molecules and the staller census (where
# native wins outright: valid ensemble in seconds vs her 0/partial after minutes).
#
#   validation/bench/speedup_report.py <genref.log> [native_times.json] [ensembles_dir]
import json
import re
import statistics
import sys
import pathlib

REPO = pathlib.Path(__file__).resolve().parent.parent.parent


def her_times(logpath, ensdir):
    """her per-molecule seconds (clean) and the staller names. Prefers the structured `_times.json`
    gen_reference now writes (robust); falls back to regex-parsing its stdout log when that's absent
    (older reference runs)."""
    times_json = ensdir / "_times.json"
    if times_json.exists():
        her, stalled = {}, set()
        for name, r in json.load(open(times_json)).items():
            if r.get("stalled"):
                stalled.add(name)
            else:
                her[name] = r["seconds"]
        return her, stalled
    return parse_genref(logpath)


def parse_genref(logpath):
    """Fallback: her per-molecule seconds (clean) and the staller names, from the generation log."""
    her, stalled = {}, set()
    # "  [  N/100] NAME  NNN confs (n=XX)   T.Ts"  and  "  [ N/100] NAME STALLED (killed ...)"
    clean_re = re.compile(r"^\s+\[\s*\d+/\d+\]\s+(\S+)\s+(\d+)\s+confs\s+\(n=\d+\)\s+([\d.]+)s")
    stall_re = re.compile(r"^\s+\[\s*\d+/\d+\]\s+(\S+)\s+STALLED")
    for line in open(logpath):
        m = clean_re.match(line)
        if m:
            her[m.group(1)] = float(m.group(3))
            continue
        m = stall_re.match(line)
        if m:
            stalled.add(m.group(1))
    return her, stalled


def main():
    if len(sys.argv) < 2:
        print("usage: speedup_report.py <genref.log> [native_times.json] [ensembles_dir]", file=sys.stderr)
        return 2
    logpath = sys.argv[1]
    ntpath = sys.argv[2] if len(sys.argv) > 2 else str(REPO / "validation/bench/native_times.json")
    ensdir = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else REPO / "validation/fixtures/divya/ensembles"

    native = {r["name"]: r for r in json.load(open(ntpath))}
    her_time, stalled_log = her_times(logpath, ensdir)
    stalls_json = ensdir / "_stalls.json"
    stalled = set(s["name"] for s in json.load(open(stalls_json))) if stalls_json.exists() else stalled_log

    # --- speedup on clean molecules (both produced an ensemble) ---
    speedups = []
    for name, ht in sorted(her_time.items()):
        nt = native.get(name, {}).get("seconds")
        if nt and nt > 0:
            speedups.append((name, ht / nt, ht, nt))
    speedups.sort(key=lambda x: x[1])
    print("=== speedup on clean molecules (her_time / native_time) ===")
    if speedups:
        s = [x[1] for x in speedups]
        med = statistics.median(s)
        print("  molecules: %d   median speedup: %.0fx   range: %.0fx .. %.0fx"
              % (len(s), med, s[0], s[-1]))
        print("  slowest-native pair: %s  her %.1fs / native %.2fs = %.0fx"
              % (speedups[0][0], speedups[0][2], speedups[0][3], speedups[0][1]))

    # --- staller census (goal dim 3): her degrades, native handles ---
    print("\n=== staller census (her pipeline stalls/fails; native's handling) ===")
    print("  molecules her pipeline stalls/fails on: %d / %d corpus" % (len(stalled), len(native)))
    nat_ok = [n for n in stalled if native.get(n, {}).get("n_conformers", 0) >= 100]
    print("  native produces a full (>=100-conf) ensemble on %d / %d of those" % (len(nat_ok), len(stalled)))
    for n in sorted(stalled):
        r = native.get(n, {})
        print("    %-14s native: %.2fs, %d confs   (her: 0/partial, >150s)"
              % (n, r.get("seconds", float("nan")), r.get("n_conformers", 0)))


if __name__ == "__main__":
    sys.exit(main())
