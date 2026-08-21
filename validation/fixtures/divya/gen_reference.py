#!/usr/bin/env python3
# Reference conformer-ensemble oracle for the Phase-3 distributional gate (CONFORMANCE.md, target B).
#
# A FAITHFUL, VERBATIM replica of Divya's ensemble block in build_ligands.py (the seeded two-stage
# core-pin recipe: seeded core embed 210185+j -> pin largest ring + exocyclic =O + exocyclic C on
# ring N -> seeded sidechain EmbedMultipleConfs). Her generation is deterministic given the SMILES
# (all seeds fixed), so this reproduces her exact ensemble. Run INSIDE her production container so it
# uses her actual RDKit build (2026.09.1pre + her amide patch):
#
#   apptainer exec --bind <repo>:<repo> <build_macrocycle_final.sif> \
#       python3 gen_reference.py <corpus.smi> <out_dir>
#
# The embed logic is copied from her code, not reimagined; only the surrounding pipeline
# (protonation/corina/db2) and the I/O are replaced. Do not "improve" the embed; fidelity is the point.
#
# STALL HANDLING: for some (molecule, core-seed) pairs the pinned coordMap is geometrically infeasible,
# and RDKit's constrained EmbedMultipleConfs thrashes maxIterations(10n)*nConf doomed attempts (~33 min)
# returning 0 conformers. Her code has no per-embed guard (only a coarse job-level scheduler signal),
# so her pipeline degrades on these too. We run each MOLECULE in a subprocess with a wall-clock budget
# and TERMINATE it if it exceeds -- a molecule that hits any stall cannot finish in the budget, so it is
# killed and recorded in _stalls.json (the reimpl-value list: where her pipeline degrades). Clean
# molecules (~30-60s) finish and get a full ensemble. This is not "her exact output" for the stallers
# (she has no clean output for them either); it is the tractable, honest handling of her defect.
import json
import multiprocessing as mp
import os
import sys
import time
from rdkit import Chem
from rdkit.Chem import AllChem

MOLECULE_TIMEOUT_S = 150  # clean molecules finish in ~30-60s; any real stall blows far past this
# Seed base for the whole recipe (core seeds = SEED_BASE+j, sidechain seed = SEED_BASE). Her production
# base is 210185; a SECOND base (env SEED_BASE=...) generates an independent ensemble whose spread vs
# the 210185 one IS the derived acceptance band for the distributional gate (her own between-seed
# variability). Default reproduces her exact output.
SEED_BASE = int(os.environ.get("SEED_BASE", "210185"))


def count_exo_rotatable(mol):  # verbatim from build_ligands.py
    ring_atoms = set(idx for ring in mol.GetRingInfo().AtomRings() for idx in ring)
    exo = 0
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        if bond.GetBondTypeAsDouble() == 1.0:
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if i in ring_atoms or j in ring_atoms:
                exo += 1
    return exo


def ensemble(smiles, name):
    """Her seeded two-stage recipe; returns a list of single-conformer Mols (rdkit_obj)."""
    mol_init = Chem.MolFromSmiles(smiles)
    if mol_init is None:
        return None
    mol_elem = Chem.AddHs(mol_init)

    exo = count_exo_rotatable(mol_init)
    if exo <= 2:
        rdkit_confs_1, sidechain_confs = 20, 10
    else:
        rdkit_confs_1, sidechain_confs = 10, 20

    sssr = Chem.GetSymmSSSR(mol_elem)
    core_atoms = list(set(max(sssr, key=len)))

    rdkit_obj = []
    for j in range(rdkit_confs_1):
        mol_init = Chem.MolFromSmiles(smiles)
        mol_elem = Chem.AddHs(mol_init)
        mol_elem.SetProp("_Name", name)
        params0 = AllChem.ETKDGv3()
        params0.useMacrocycleTorsions = True
        params0.useExpTorsionAnglePrefs = True
        params0.randomSeed = SEED_BASE + j
        result = AllChem.EmbedMolecule(mol_elem, params0)
        if result == -1:
            continue
        cmap = {}
        for k in core_atoms:
            cmap[k] = mol_elem.GetConformer().GetAtomPosition(k)
            atom = mol_elem.GetAtomWithIdx(k)
            for neighbor in atom.GetNeighbors():
                nb_idx = neighbor.GetIdx()
                bond = mol_elem.GetBondBetweenAtoms(k, nb_idx)
                if (neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0 and nb_idx not in cmap):
                    cmap[nb_idx] = mol_elem.GetConformer().GetAtomPosition(nb_idx)
                if (atom.GetAtomicNum() == 7 and neighbor.GetAtomicNum() == 6 and nb_idx not in core_atoms and nb_idx not in cmap):
                    cmap[nb_idx] = mol_elem.GetConformer().GetAtomPosition(nb_idx)
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = SEED_BASE
        params.useExpTorsionAnglePrefs = True
        params.SetCoordMap(cmap)
        AllChem.EmbedMultipleConfs(mol_elem, sidechain_confs, params)
        for l in range(mol_elem.GetNumConformers()):
            conf = mol_elem.GetConformer(l)
            conf_mol = Chem.Mol(mol_elem)
            conf_mol.RemoveAllConformers()
            conf_mol.AddConformer(conf, assignId=True)
            rdkit_obj.append(conf_mol)
    return rdkit_obj


def _worker(smiles, name, out_path):
    """Runs one molecule's full ensemble and writes the .confs.json. Killed by the parent on timeout."""
    confs = ensemble(smiles, name)
    if not confs:
        return
    n_atoms = confs[0].GetNumAtoms()
    # flattened coords per conformer, RDKit AddHs order (= native spec atom order; heavy indices align)
    flat = [[float(v) for row in m.GetConformer().GetPositions() for v in row] for m in confs]
    rec = {"name": name, "n_conformers": len(confs), "n_atoms": n_atoms, "conformers": flat}
    with open(out_path, "w") as w:
        json.dump(rec, w)


def main():
    if len(sys.argv) != 3:
        print("usage: gen_reference.py <corpus.smi> <out_dir>", file=sys.stderr)
        return 2
    corpus, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    nproc = int(os.environ.get("GENREF_NPROC", "24"))
    print("rdkit", Chem.rdBase.rdkitVersion, "molecule_timeout=%ds nproc=%d" % (MOLECULE_TIMEOUT_S, nproc))
    ctx = mp.get_context("fork")

    mols = []
    for i, line in enumerate(open(corpus)):
        tok = line.split()
        if tok:
            mols.append((tok[0], tok[1] if len(tok) > 1 else "mol%04d" % i))

    # Bounded-concurrency pool: up to `nproc` molecules run at once, each in its own timed subprocess
    # (per-molecule wall-clock budget, terminated on stall). Molecules are independent + deterministic,
    # so concurrency changes only throughput, not any molecule's output.
    n_ok, stalls, done = 0, [], 0
    times = {}    # name -> {seconds, confs, n_atoms, stalled} — the structured form of the human log,
                  # so consumers (speedup_report.py) read a JSON instead of regex-parsing stdout.
    active = {}   # proc -> (name, smiles, out_path, t0)
    idx = 0
    while idx < len(mols) or active:
        while idx < len(mols) and len(active) < nproc:
            smiles, name = mols[idx]
            idx += 1
            out_path = os.path.join(out_dir, "%s.confs.json" % name)
            if os.path.exists(out_path):
                os.remove(out_path)
            p = ctx.Process(target=_worker, args=(smiles, name, out_path))
            p.start()
            active[p] = (name, smiles, out_path, time.time())
        for p in list(active):
            name, smiles, out_path, t0 = active[p]
            if not p.is_alive():
                p.join()
                done += 1
                secs = time.time() - t0
                if os.path.exists(out_path):
                    d = json.load(open(out_path))
                    n_ok += 1
                    times[name] = {"seconds": secs, "confs": d["n_conformers"],
                                   "n_atoms": d["n_atoms"], "stalled": False}
                    print("  [%3d/%d] %-14s %3d confs (n=%d)  %5.1fs"
                          % (done, len(mols), name, d["n_conformers"], d["n_atoms"], secs))
                else:
                    stalls.append({"name": name, "smiles": smiles, "reason": "parse/embed fail"})
                    times[name] = {"seconds": secs, "confs": 0, "n_atoms": 0, "stalled": True}
                    print("  [%3d/%d] %-14s SKIP (parse/embed fail)" % (done, len(mols), name))
                del active[p]
            elif time.time() - t0 > MOLECULE_TIMEOUT_S:
                p.terminate()
                p.join()
                done += 1
                stalls.append({"name": name, "smiles": smiles, "reason": "timeout>%ds" % MOLECULE_TIMEOUT_S})
                times[name] = {"seconds": time.time() - t0, "confs": 0, "n_atoms": 0, "stalled": True}
                print("  [%3d/%d] %-14s STALLED (killed >%ds) - her pipeline degrades here"
                      % (done, len(mols), name, MOLECULE_TIMEOUT_S))
                del active[p]
        sys.stdout.flush()
        time.sleep(0.3)
    json.dump(stalls, open(os.path.join(out_dir, "_stalls.json"), "w"), indent=1)
    json.dump(times, open(os.path.join(out_dir, "_times.json"), "w"), indent=1)
    print("wrote %d clean ensembles; %d molecules her pipeline stalls/fails on (see _stalls.json)"
          % (n_ok, len(stalls)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
