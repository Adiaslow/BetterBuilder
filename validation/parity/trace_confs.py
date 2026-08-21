"""Capture the oracle's conformer ensembles by tracing its unmodified execution.

build_ligands.py writes its conformers to a temporary SDF and deletes it, so the ensemble is never
an output of hers we could simply keep. Rather than intercept that temp file — or reimplement the
recipe — this runs her file unmodified under sys.settrace and reads `rdkit_obj` out of the live
frame once the seed loop has filled it.

    python3 trace_confs.py <protomers.smi> <outdir>

<protomers.smi> is `<smiles> <name>`; it is fed to her file on stdin, which is how the production
build_ligands.py reads protomers. /data/input.smi must also be bound (her neutral-SMILES dict).

Writes <outdir>/<molecule>.confs.json: one entry per conformer, coordinates in RDKit atom order.
"""
import json
import os
import runpy
import sys
import traceback

BUILD = "/soft/DOCK-macro-latest6/ligand/generate/build_ligands.py"
# The statement after which rdkit_obj holds the full ensemble. Located by source text rather than
# by a fixed line number, which differs between builds of her file.
CAPTURE_MARKER = "tmp_out_sdf = tempfile.NamedTemporaryFile"


def _capture_line(path, marker):
    with open(path) as f:
        for i, line in enumerate(f, 1):
            if marker in line:
                return i
    raise SystemExit(f"trace_confs: marker {marker!r} not found in {path}")


CAPTURE_LINE = _capture_line(BUILD, CAPTURE_MARKER)

_seen = {}


_BUILD_NAMES = {BUILD, os.path.realpath(BUILD), os.path.basename(BUILD)}
_seen_files = set()


def _is_build(frame):
    co = frame.f_code.co_filename
    return co in _BUILD_NAMES or os.path.basename(co) == os.path.basename(BUILD)


def tracer(frame, event, arg):
    # Only trace inside her file. Returning None for other frames stops line-tracing from
    # descending into RDKit, which otherwise dominates runtime by orders of magnitude.
    if not _is_build(frame):
        if os.environ.get("BB_TRACE_DEBUG") and event == "call":
            _seen_files.add(frame.f_code.co_filename)
        return None
    if event == "line" and frame.f_lineno == CAPTURE_LINE:
        obj = frame.f_locals.get("rdkit_obj")
        mol = frame.f_locals.get("mol")
        name = getattr(mol, "name", None)
        if obj and name and name not in _seen:
            confs = []
            for cm in obj:
                c = cm.GetConformer()
                confs.append([v for i in range(cm.GetNumAtoms())
                              for v in (c.GetAtomPosition(i).x,
                                        c.GetAtomPosition(i).y,
                                        c.GetAtomPosition(i).z)])
            _seen[name] = {
                "name": name,
                "smiles": getattr(mol, "smiles", None),
                "n_atoms": obj[0].GetNumAtoms(),
                "n_conformers": len(confs),
                "conformers": confs,
            }
            print(f"  {name}: {len(confs)} conformers x {obj[0].GetNumAtoms()} atoms")
    return tracer


def main(protomers, outdir):
    # The production build_ligands.py reads protomers from stdin as `<smiles> <name>` and assigns
    # protomer ids itself; it also loads /data/input.smi as a neutral-SMILES dictionary. (The
    # commented-out sys.argv path in that file is the older contract.)
    sys.path.insert(0, os.path.dirname(BUILD))
    sys.argv = [BUILD]
    sys.stdin = open(protomers)
    os.makedirs(outdir, exist_ok=True)
    sys.settrace(tracer)
    try:
        runpy.run_path(BUILD, run_name="__main__")
    except BaseException as e:
        traceback.print_exc()
        print(f"  (pipeline ended: {type(e).__name__})")
    finally:
        sys.settrace(None)
    for name, data in _seen.items():
        with open(os.path.join(outdir, f"{name.split('_')[0]}.confs.json"), "w") as f:
            json.dump(data, f)
    print(f"  wrote {len(_seen)} ensembles to {outdir}")
    if not _seen:
        print(f"  capture line {CAPTURE_LINE} ({CAPTURE_MARKER!r}) never matched in {BUILD}")
        if os.environ.get("BB_TRACE_DEBUG"):
            for f in sorted(_seen_files)[:20]:
                print(f"    saw frame from {f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
