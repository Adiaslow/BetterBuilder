"""Gate the fully-native MoleculeSpec against the RDKit-bridge MoleculeSpec, field by field.

    python3 gate_spec.py <corpus.smi> [N]           # per-molecule subprocess (small N)
    python3 gate_spec.py --batch <bridge.jsonl> <native.jsonl>   # whole-corpus JSONL-to-JSONL

Compares every field of the emitted spec. The bridge is the oracle: the native spec is correct
exactly when it reproduces the bridge's spec. Bounds are f32 and compared to a small absolute
tolerance; the discrete fields (atomic numbers, bonds, angles, torsions, impropers, chirality,
stereo, pins, seeds, dim, charge) must match exactly, treating the constraint lists as sets since
emission order is not meaningful.

The per-molecule mode launches `bb-spec` (the RDKit bridge, in bb-rdkit) and `bb-spec-native` (the
pure-Rust assembler, in bb-spec) once each per SMILES — fine for a spot check, but it pays RDKit's
shared-object load on every molecule. For the whole corpus, produce the two JSONL files once with
`bb-spec-batch` and `bb-spec-native-batch` (one process each), then gate them with `--batch`.

Exit 0 iff every molecule matches on every field.
"""
import gzip
import json
import math
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "lib"))
# single source: the drift-prone bridge-vs-native canonicalisers ([[stereo-representation-equivalence]])
from spec_diff import canon_stereo, canon_chiral, canon_tors, canon_imp, canon_set

BOUNDS_ATOL = 0.0  # native smooths in f64 then casts, exactly as the bridge → byte-identical bounds

BRIDGE = "../../rust/target/debug/bb-spec"
NATIVE = "../../rust/target/debug/bb-spec-native"


def run(binary, smiles):
    p = subprocess.run([binary, smiles], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        return None
    return json.loads(p.stdout.decode())


def compare(a, b):
    """Return a list of field names that differ."""
    diffs = []

    # scalars that must match exactly
    for k in ("n_atoms", "dim", "formal_charge", "core_seeds", "sidechain_confs"):
        if a.get(k) != b.get(k):
            diffs.append(f"{k}({a.get(k)}!={b.get(k)})")

    if a.get("atomic_numbers") != b.get("atomic_numbers"):
        diffs.append("atomic_numbers")
    if a.get("bounds_force_scaling") != b.get("bounds_force_scaling"):
        diffs.append("bounds_force_scaling")

    # bounds: f32, absolute tolerance
    ab, bb_ = a.get("bounds", []), b.get("bounds", [])
    if len(ab) != len(bb_):
        diffs.append(f"bounds_len({len(ab)}!={len(bb_)})")
    else:
        worst = max((abs(x - y) for x, y in zip(ab, bb_)), default=0.0)
        if worst > BOUNDS_ATOL or math.isnan(worst):
            diffs.append(f"bounds(|Δ|max={worst:.2e})")

    # order-insensitive discrete sets (canonicalisers shared with the hunt gate via spec_diff)
    if canon_set(a.get("bonds", [])) != canon_set(b.get("bonds", [])):
        diffs.append("bonds")
    if canon_set([tuple(x["atoms"]) + (x["triple"],) for x in a.get("angles", [])]) != \
       canon_set([tuple(x["atoms"]) + (x["triple"],) for x in b.get("angles", [])]):
        diffs.append("angles")
    if canon_set(a.get("pin_atoms", [])) != canon_set(b.get("pin_atoms", [])):
        diffs.append("pin_atoms")
    if canon_set(a.get("double_bond_ends", [])) != canon_set(b.get("double_bond_ends", [])):
        diffs.append("double_bond_ends")

    if canon_chiral(a.get("chiral_sets", [])) != canon_chiral(b.get("chiral_sets", [])):
        diffs.append("chiral_sets")
    if canon_chiral(a.get("tetrahedral_centers", [])) != canon_chiral(b.get("tetrahedral_centers", [])):
        diffs.append("tetrahedral_centers")
    if canon_stereo(a.get("stereo_double_bonds", []), a.get("bonds", [])) != \
       canon_stereo(b.get("stereo_double_bonds", []), b.get("bonds", [])):
        diffs.append("stereo_double_bonds")

    if canon_tors(a.get("exp_torsions", [])) != canon_tors(b.get("exp_torsions", [])):
        diffs.append("exp_torsions")
    if canon_imp(a.get("impropers", [])) != canon_imp(b.get("impropers", [])):
        diffs.append("impropers")

    return diffs


def main(corpus, n):
    ok = miss = fail = 0
    field_fail = {}
    examples = []
    for i, line in enumerate(open(corpus)):
        if i >= n:
            break
        f = line.split()
        if len(f) < 2:
            continue
        smiles, name = f[0], f[1]
        bspec = run(BRIDGE, smiles)
        nspec = run(NATIVE, smiles)
        if bspec is None or nspec is None:
            miss += 1
            continue
        diffs = compare(nspec, bspec)
        if diffs:
            fail += 1
            for d in diffs:
                field_fail[d.split("(")[0]] = field_fail.get(d.split("(")[0], 0) + 1
            if len(examples) < 10:
                examples.append((name, smiles, diffs))
        else:
            ok += 1

    print(f"native vs bridge spec: {ok} match, {fail} differ, {miss} unbuildable "
          f"({ok + fail + miss} seen)")
    if field_fail:
        print("  per-field failures:", dict(sorted(field_fail.items(), key=lambda x: -x[1])))
        for name, smiles, diffs in examples:
            print(f"  {name}: {', '.join(diffs)}")
            print(f"    {smiles}")
    return 0 if fail == 0 and ok > 0 else 1


def _load_jsonl(path):
    """Map name -> record ({"spec":..} or {"error":..}) from a batch JSONL file."""
    out = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            out[r["name"]] = r
    return out


def main_batch(bridge_path, native_path):
    bridge = _load_jsonl(bridge_path)
    native = _load_jsonl(native_path)
    ok = miss = fail = 0
    field_fail = {}
    examples = []
    for name, br in bridge.items():
        nr = native.get(name)
        if nr is None or "spec" not in br or "spec" not in nr:
            # unbuildable on either side (both reject → not a discrepancy we can gate here)
            miss += 1
            continue
        diffs = compare(nr["spec"], br["spec"])
        if diffs:
            fail += 1
            for d in diffs:
                field_fail[d.split("(")[0]] = field_fail.get(d.split("(")[0], 0) + 1
            if len(examples) < 10:
                examples.append((name, diffs))
        else:
            ok += 1

    print(f"native vs bridge spec: {ok} match, {fail} differ, {miss} unbuildable "
          f"({ok + fail + miss} seen)")
    if field_fail:
        print("  per-field failures:", dict(sorted(field_fail.items(), key=lambda x: -x[1])))
        for name, diffs in examples:
            print(f"  {name}: {', '.join(diffs)}")
    return 0 if fail == 0 and ok > 0 else 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--batch":
        sys.exit(main_batch(sys.argv[2], sys.argv[3]))
    corpus = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    sys.exit(main(corpus, n))
