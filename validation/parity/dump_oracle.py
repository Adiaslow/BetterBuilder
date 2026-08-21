"""Dump stage-level state from the oracle's own modules, for differential testing.

Runs inside the pipeline container and imports Divya's modules **unmodified**, calling them as a
library. Nothing here reimplements or edits her code; every value written out is one her code
computed. BetterBuilder's Rust is then checked against these values stage by stage, so a divergence
identifies which stage broke rather than only that the final artifact differs.

    apptainer exec -B ... <image> python3 dump_oracle.py <mol2> <solv> <out.json>

The mol2 and solv are the oracle's own intermediates, taken from a pipeline run.
"""
import hashlib
import io
import json
import sys

sys.path[:0] = [
    "/soft/DOCK-macro-latest6/ligand/mol2db2_py3_strain",
    "/soft/DOCK-macro-latest6/ligand/strain",
]

import clash  # noqa: E402
import hierarchy  # noqa: E402
import mol2  # noqa: E402
import solv  # noqa: E402

CLASHFILE = "/soft/DOCK-macro-latest6/ligand/mol2db2/clashfile.txt"
TOLERANCE = 0.001  # mol2db2.py's --disttol default


def jsonable(v):
    """Reduce a value to something JSON can hold, preserving what we compare on."""
    if isinstance(v, (int, float, str, bool)) or v is None:
        return v
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    if isinstance(v, (set, frozenset)):
        return sorted(jsonable(x) for x in v)
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    return repr(v)


def main(mol2path, solvpath, outpath):
    mol2data = mol2.Mol2(mol2path)
    mol2data.convertDockTypes(None)  # None = the built-in table, as mol2db2_quick uses
    mol2data.addColors(None)
    solvdata = solv.Solv(solvpath)
    clashDecider = clash.Clash(CLASHFILE)

    out = {
        "source": {"mol2": mol2path, "solv": solvpath, "tolerance": TOLERANCE},
        "mol2": {
            "name": mol2data.name,
            "n_atoms": len(mol2data.atomName),
            "n_conformers": len(mol2data.atomXyz),
            "atom_name": jsonable(mol2data.atomName),
            "atom_type": jsonable(mol2data.atomType),
            "dock_num": jsonable(getattr(mol2data, "dockNum", None)),
            "color_num": jsonable(getattr(mol2data, "colorNum", None)),
            "atom_bonds": jsonable(mol2data.atomBonds),
        },
        "solv": {
            "charge": jsonable(solvdata.charge),
            "polar_solv": jsonable(solvdata.polarSolv),
            "surface": jsonable(solvdata.surface),
            "apolar_solv": jsonable(solvdata.apolarSolv),
            "solv": jsonable(solvdata.solv),
        },
    }

    h = hierarchy.Hierarchy(mol2data, clashDecider, tolerance=TOLERANCE, solvdata=solvdata)

    # Every non-callable, non-private attribute the hierarchy exposes after construction. Captured
    # wholesale rather than hand-picked, so a stage we have not thought about is still covered.
    hier = {}
    for attr in sorted(dir(h)):
        if attr.startswith("__"):
            continue
        v = getattr(h, attr)
        if callable(v) or attr in ("mol2data", "solvdata"):
            continue
        hier[attr] = jsonable(v)
    out["hierarchy"] = hier

    sio = io.StringIO()
    h.writeFile(fileHandle=sio)
    db2 = sio.getvalue()
    out["db2"] = {
        "bytes": len(db2),
        "sha256": hashlib.sha256(db2.encode()).hexdigest(),
        "first_lines": db2.split("\n")[:12],
    }

    with open(outpath, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(f"{mol2data.name}: {out['mol2']['n_atoms']} atoms, "
          f"{out['mol2']['n_conformers']} confs, db2 {out['db2']['bytes']} bytes")


if __name__ == "__main__":
    main(*sys.argv[1:4])
