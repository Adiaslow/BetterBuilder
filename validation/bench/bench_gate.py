#!/usr/bin/env python3
# Failable speed gate over validation/bench/run.sh output. Encodes the two invariants the bench exists
# to protect (see [[speed-floor-is-cpp-rdkit]]):
#
#   1. EQUAL WORK — the Rust and C++ harness for a component must report the SAME work quantity (atom /
#      torsion / constrained-pair counts exact; content checksums within a small relative tolerance).
#      Without this the ms/mol are not comparable and any speed claim is meaningless.
#   2. C++ IS THE FLOOR — the pure-Rust component must be at least as fast as C++ RDKit doing that same
#      work (that is the whole justification for replacing the RDKit routine).
#
# Reads run.sh's stdout on stdin (or a file arg). Exits non-zero if any component violates either
# invariant, or if a component is missing its Rust or C++ half. Deterministic on the equal-work check;
# the speed check compares the harnesses' own warmed-up aggregate timings.
import re
import sys

# Rust must be under this multiple of the C++ time. 1.0 = "must beat the floor"; the components clear it
# by 1.2-3.0x today, so a strict floor won't flake on the aggregate-of-100+ timings the harnesses report.
SPEED_MARGIN = float(__import__("os").environ.get("BB_BENCH_SPEED_MARGIN", "1.0"))
CHECKSUM_RTOL = 1e-5  # Rust prints full f64, C++ prints fewer digits of the same sum

MSPER = re.compile(r"=\s*([\d.]+)\s*ms/mol")
COUNTS = re.compile(r"(atoms_typed|torsions|constrained_pairs)\s+(\d+)")
CHECKSUM = re.compile(r"checksum\s+([\d.]+)")


def parse(text):
    """{component: {'rust': (ms, counts, checksum), 'cpp': (...)}} from run.sh output."""
    out, comp, side = {}, None, None
    for line in text.splitlines():
        s = line.strip()
        m = re.match(r"===\s*(.+?)\s*===", s)
        if m:
            comp = m.group(1)
            out.setdefault(comp, {})
            continue
        if s.startswith("--- Rust"):
            side = "rust"
            continue
        if s.startswith("--- C++") or s.startswith("--- Cpp"):
            side = "cpp"
            continue
        mm = MSPER.search(line)
        if mm and comp and side:
            ms = float(mm.group(1))
            counts = {k: int(v) for k, v in COUNTS.findall(line)}
            cs = CHECKSUM.search(line)
            out[comp][side] = (ms, counts, float(cs.group(1)) if cs else None)
    return out


def main():
    text = open(sys.argv[1]).read() if len(sys.argv) > 1 else sys.stdin.read()
    comps = parse(text)
    if not comps:
        print("bench gate: no components parsed from run.sh output", file=sys.stderr)
        return 2
    failed = []
    print("=== bench gate: Rust component vs C++ RDKit floor (equal work) ===")
    for comp, sides in sorted(comps.items()):
        if "rust" not in sides or "cpp" not in sides:
            print(f"  {comp:22s} MISSING {'rust' if 'rust' not in sides else 'cpp'} half")
            failed.append(comp)
            continue
        (r_ms, r_ct, r_cs), (c_ms, c_ct, c_cs) = sides["rust"], sides["cpp"]
        # invariant 1: equal work
        work_ok = r_ct == c_ct
        if r_cs is not None and c_cs is not None:
            denom = max(abs(r_cs), abs(c_cs), 1.0)
            work_ok = work_ok and abs(r_cs - c_cs) / denom <= CHECKSUM_RTOL
        # invariant 2: Rust beats the C++ floor
        speed_ok = r_ms <= c_ms * SPEED_MARGIN
        ratio = c_ms / r_ms if r_ms > 0 else float("inf")
        tag = "OK" if (work_ok and speed_ok) else "FAIL"
        if not work_ok:
            tag += " (work mismatch: rust %s/%s vs cpp %s/%s)" % (r_ct, r_cs, c_ct, c_cs)
        if not speed_ok:
            tag += " (rust %.4f > floor %.4f ms/mol)" % (r_ms, c_ms)
        print(f"  {comp:22s} rust {r_ms:8.4f}  cpp {c_ms:8.4f} ms/mol  = {ratio:4.2f}x  {tag}")
        if not (work_ok and speed_ok):
            failed.append(comp)
    if failed:
        print(f"  => BENCH GATE FAILED on {len(failed)}: {', '.join(failed)}", file=sys.stderr)
        return 1
    print(f"  => all {len(comps)} components beat the C++ floor doing equal work")
    return 0


if __name__ == "__main__":
    sys.exit(main())
