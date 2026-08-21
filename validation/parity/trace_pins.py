"""Capture the oracle's coordMap (pinned atoms) by tracing its unmodified execution.

`cmap` is built inline inside build_ligands.py's seed loop, so unlike count_exo_rotatable it has no
extraction boundary. Rather than reimplement it, this runs her file unmodified under sys.settrace
and reads `cmap` out of the live frame at the line where she hands it to SetCoordMap. Every value
recorded is one her code computed.

The same frame also holds the recipe parameters, so pins and `rdkit_confs_1`/`sidechain_confs` come
from one traced execution rather than two capture paths.

    python3 trace_pins.py <protomers.smi> <out.json>

<protomers.smi> is `<smiles> <name>`; it is fed to her file on stdin, which is how the production
build_ligands.py reads protomers. /data/input.smi must also be bound (her neutral-SMILES dict).
"""
import json
import os
import runpy
import sys

BUILD = "/soft/DOCK-macro-latest6/ligand/generate/build_ligands.py"
# Located by source text rather than a fixed line number, which differs between builds of her file.
CAPTURE_MARKER = "params.SetCoordMap(cmap)"


def _capture_line(path, marker):
    with open(path) as f:
        for i, line in enumerate(f, 1):
            if marker in line:
                return i
    raise SystemExit(f"trace_pins: marker {marker!r} not found in {path}")


CAPTURE_LINE = _capture_line(BUILD, CAPTURE_MARKER)
_BUILD_NAMES = {BUILD, os.path.realpath(BUILD), os.path.basename(BUILD)}

captured = []


def _is_build(frame):
    co = frame.f_code.co_filename
    return co in _BUILD_NAMES or os.path.basename(co) == os.path.basename(BUILD)


def tracer(frame, event, arg):
    # Only trace inside her file. Returning None for other frames stops line-tracing from
    # descending into RDKit, which otherwise dominates runtime by orders of magnitude.
    if not _is_build(frame):
        return None
    if event == "line" and frame.f_lineno == CAPTURE_LINE:
        cmap = frame.f_locals.get("cmap")
        mol = frame.f_locals.get("mol")
        if cmap is not None:
            captured.append({
                "name": getattr(mol, "name", None),
                "pins": sorted(int(k) for k in cmap.keys()),
                "core_atoms": sorted(int(x) for x in frame.f_locals.get("core_atoms", [])),
                "rdkit_confs_1": frame.f_locals.get("rdkit_confs_1"),
                "sidechain_confs": frame.f_locals.get("sidechain_confs"),
            })
    return tracer


def main(protomers, outpath):
    # `python3 build_ligands.py` puts the script's directory on sys.path; runpy does not.
    sys.path.insert(0, os.path.dirname(BUILD))
    sys.argv = [BUILD]
    sys.stdin = open(protomers)
    sys.settrace(tracer)
    try:
        runpy.run_path(BUILD, run_name="__main__")
    except BaseException as e:            # the run ends at shutil.move to /data; we have what we need
        print(f"  (pipeline ended: {type(e).__name__})")
    finally:
        sys.settrace(None)
    # one entry per core seed; the pin set is identical across seeds for a molecule
    by_mol = {}
    for c in captured:
        by_mol.setdefault(c["name"], []).append(c)
    out = {}
    for name, entries in by_mol.items():
        pins = {tuple(e["pins"]) for e in entries}
        out[name] = {
            "n_seeds_traced": len(entries),
            "pin_sets_distinct": len(pins),
            "pins": sorted(entries[0]["pins"]),
            "core_atoms": entries[0]["core_atoms"],
            "rdkit_confs_1": entries[0]["rdkit_confs_1"],
            "sidechain_confs": entries[0]["sidechain_confs"],
        }
        print(f"  {name}: {len(entries[0]['pins'])} pinned atoms, "
              f"{len(entries[0]['core_atoms'])} core, {len(entries)} seeds traced")
    with open(outpath, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
