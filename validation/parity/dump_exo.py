"""Run her `count_exo_rotatable` over a corpus by extracting it from her file.

The function is nested inside build_ligands.py's per-molecule loop, so it cannot be imported. This
locates its `ast.FunctionDef`, compiles that node alone and executes it — her exact bytes, with no
reimplementation and without running the pipeline. It needs no embedding, so it covers thousands of
molecules in seconds rather than the ~1 hour per 25 that tracing costs.

    python3 dump_exo.py <corpus.smi> <out.json>

Writes {name: {"exo": n, "core_seeds": s, "sidechain_confs": c}} using her branch thresholds, read
from her source rather than restated here.
"""
import ast
import json
import re
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BUILD = "/soft/DOCK-macro-latest6/ligand/generate/build_ligands.py"
FUNC = "count_exo_rotatable"


def extract(path, name):
    """Compile and execute just the named FunctionDef from her file."""
    src = open(path).read()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            mod = ast.Module(body=[node], type_ignores=[])
            ast.fix_missing_locations(mod)
            ns = {"Chem": Chem}
            exec(compile(mod, path, "exec"), ns)  # noqa: S102 — her bytes, deliberately
            return ns[name]
    raise SystemExit(f"dump_exo: {name} not found in {path}")


def branch_thresholds(path):
    """Read the exo -> (core_seeds, sidechain_confs) mapping out of her source."""
    src = open(path).read()
    m = re.search(
        r"if exo\s*<=\s*(\d+):\s*\n\s*rdkit_confs_1\s*=\s*(\d+)\s*\n\s*sidechain_confs\s*=\s*(\d+)"
        r"\s*\n\s*elif exo\s*>=\s*(\d+):\s*\n\s*rdkit_confs_1\s*=\s*(\d+)\s*\n\s*sidechain_confs\s*=\s*(\d+)",
        src,
    )
    if not m:
        raise SystemExit("dump_exo: could not read the branch thresholds from her source")
    lo, a1, b1, hi, a2, b2 = (int(x) for x in m.groups())
    return lo, (a1, b1), hi, (a2, b2)


def main(corpus, outpath):
    count_exo = extract(BUILD, FUNC)
    lo, low_branch, hi, high_branch = branch_thresholds(BUILD)
    out, bad = {}, 0
    for line in open(corpus):
        f = line.split()
        if len(f) < 2:
            continue
        # she calls it on the pre-addHs molecule from SMILES
        mol = Chem.MolFromSmiles(f[0])
        if mol is None:
            bad += 1
            continue
        exo = count_exo(mol)
        seeds, side = low_branch if exo <= lo else high_branch
        out[f[1]] = {"exo": exo, "core_seeds": seeds, "sidechain_confs": side}
    with open(outpath, "w") as fh:
        json.dump(out, fh, sort_keys=True)
    branches = {}
    for v in out.values():
        k = (v["core_seeds"], v["sidechain_confs"])
        branches[k] = branches.get(k, 0) + 1
    print(f"  {len(out)} molecules ({bad} unparseable); thresholds exo<={lo} -> {low_branch}, "
          f"exo>={hi} -> {high_branch}")
    for k, n in sorted(branches.items()):
        print(f"    {k[0]}x{k[1]}: {n}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
