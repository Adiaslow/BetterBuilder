#!/usr/bin/env python3
# Divergence-hunting differential (goal criterion 5). For every generated molecule, run the RDKit
# BRIDGE (bb-spec, the oracle) and NATIVE (bb-spec-native) and diff the MoleculeSpec field-by-field.
# A divergence on a HARD field (elements/AddHs order, bonds, recipe counts, torsion/improper topology)
# is a perception bug or a spec gap; a NUMERIC field is compared with the vetted tolerance; STEREO
# fields are reported separately because some differ only in REPRESENTATION ([[stereo-representation-
# equivalence]]) — a raw mismatch there is a finding to inspect, not necessarily a defect.
#
# "no diffs" only means "searched hard" if coverage is measured, so this reports how many molecules
# were hunted and the structural range, alongside any divergence. FAILABLE: exits nonzero if any HARD
# divergence is found.
#   apptainer exec --bind <repo>:<repo> <sif> python3 validation/hunt/differential.py <generated.smi>
import collections
import json
import multiprocessing as mp
import os
import subprocess
import sys
import pathlib

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "validation/lib"))
from spec_diff import canon_stereo, num_diff, tup  # single source: the drift-prone canonicalisers

BRIDGE = REPO / "rust/target/release/bb-spec"
NATIVE = REPO / "rust/target/release/bb-spec-native"
NUM_TOL = 1e-4   # f32 bounds storage; native vs bridge agree well within this


def spec(binary, smi):
    r = subprocess.run([str(binary), smi], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                       universal_newlines=True)
    if r.returncode != 0 or not r.stdout.strip():
        return None
    try:
        return json.loads(r.stdout)
    except json.JSONDecodeError:
        return None


def compare(bridge, native):
    """Return (hard_diffs, stereo_diffs, numeric_diffs) as lists of field names."""
    hard, stereo, numeric = [], [], []
    # HARD structural / recipe
    for k in ("n_atoms", "dim", "formal_charge", "core_seeds", "sidechain_confs"):
        if bridge.get(k) != native.get(k):
            hard.append(k)
    if bridge.get("atomic_numbers") != native.get("atomic_numbers"):
        hard.append("atomic_numbers")            # elements + AddHs order (index alignment)
    if set(map(frozenset, bridge.get("bonds", []))) != set(map(frozenset, native.get("bonds", []))):
        hard.append("bonds")
    if set(bridge.get("pin_atoms") or []) != set(native.get("pin_atoms") or []):
        hard.append("pin_atoms")
    # torsion / improper TOPOLOGY (hard): the set of atom-tuples must match (amide patch lives here)
    for k in ("exp_torsions", "impropers", "angles"):
        bs = set(tup(x) for x in bridge.get(k, []) or [])
        ns = set(tup(x) for x in native.get(k, []) or [])
        if bs != ns:
            hard.append(k + "_topology")
    # NUMERIC — only `bounds` (the smoothed f32 matrix) is emitted by BOTH; raw_bounds is native-only.
    if num_diff(bridge.get("bounds"), native.get("bounds"), NUM_TOL):
        numeric.append("bounds")
    # STEREO (representation-sensitive: report, don't hard-fail)
    for k in ("chiral_sets", "tetrahedral_centers", "double_bond_ends"):
        bs = set(tup(x) for x in bridge.get(k, []) or [])
        ns = set(tup(x) for x in native.get(k, []) or [])
        if bs != ns:
            stereo.append(k)
    # stereo double bonds: compare on the REPRESENTATION-INVARIANT canonical form (proves the two
    # encodings are the same geometry rather than flagging a different-reference-atom convention).
    bb = canon_stereo(bridge.get("stereo_double_bonds", []) or [], bridge.get("bonds", []) or [])
    nb = canon_stereo(native.get("stereo_double_bonds", []) or [], native.get("bonds", []) or [])
    if bb != nb:
        stereo.append("stereo_double_bonds")
    return hard, stereo, numeric


def hunt_one(smi):
    """Diff one molecule; returns (kind, hard, stereo, numeric, smi). kind: 'parse'|'both_reject'|'ok'."""
    b, nat = spec(BRIDGE, smi), spec(NATIVE, smi)
    if (b is None) != (nat is None):
        return ("parse", [f"bridge={'ok' if b else 'fail'} native={'ok' if nat else 'fail'}"], [], [], smi)
    if b is None:
        return ("both_reject", [], [], [], smi)
    hard, stereo, numeric = compare(b, nat)
    return ("ok", hard, stereo, numeric, smi)


def main():
    gen = sys.argv[1] if len(sys.argv) > 1 else str(REPO / "validation/hunt/generated.smi")
    mols = [l.split()[0] for l in open(gen) if l.split()]
    hard_ct, stereo_ct, num_ct = collections.Counter(), collections.Counter(), collections.Counter()
    parse_div, n = 0, 0
    hard_ex, stereo_ex = [], []
    nproc = int(os.environ.get("HUNT_NPROC", "24"))
    with mp.get_context("fork").Pool(nproc) as pool:
        results = pool.imap_unordered(hunt_one, mols, chunksize=8)
        for kind, hard, stereo, numeric, smi in results:
            if kind == "parse":
                parse_div += 1
                if len(hard_ex) < 10:
                    hard_ex.append(f"PARSE {hard[0]}: {smi}")
                continue
            if kind == "both_reject":
                continue
            n += 1
            for f in hard:
                hard_ct[f] += 1
            for f in stereo:
                stereo_ct[f] += 1
            for f in numeric:
                num_ct[f] += 1
            if hard and len(hard_ex) < 10:
                hard_ex.append(f"HARD {hard}: {smi}")
            if stereo and not hard and len(stereo_ex) < 6:
                stereo_ex.append(f"STEREO {stereo}: {smi}")

    print(f"=== divergence hunt: {n} molecules diffed (native vs RDKit bridge), {len(mols)} generated ===")
    print(f"  PARSE divergences (one accepts, other rejects): {parse_div}")
    print(f"  HARD    divergences: {dict(hard_ct) or 'NONE'}")
    print(f"  NUMERIC divergences (> {NUM_TOL:.0e}): {dict(num_ct) or 'NONE'}")
    print(f"  STEREO  divergences (representation — inspect): {dict(stereo_ct) or 'NONE'}")
    for e in hard_ex:
        print("   ", e)
    for e in stereo_ex:
        print("   ", e)
    total_hard = parse_div + sum(hard_ct.values())
    print(f"\n  {'FAIL' if total_hard else 'PASS'}: {total_hard} hard/parse divergences over {n} molecules")
    return 1 if total_hard else 0


if __name__ == "__main__":
    sys.exit(main())
