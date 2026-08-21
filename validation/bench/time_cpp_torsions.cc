// C++ baseline for torsion assignment, against the same patched RDKit the bridge links.
//
//   time_cpp_torsions <corpus.smi>
//
// The speed floor for bb-perceive's torsion stage is C++ RDKit, not the Python binding (which adds
// per-call interpreter overhead). This does end-to-end what bb-torsions does: SMILES -> AddHs ->
// getExperimentalTorsions, timed over the whole corpus after a warm-up, so the two numbers compare
// like for like. It calls ForceFields::CrystalFF::getExperimentalTorsions directly rather than the
// Python GetExperimentalTorsions wrapper.
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: time_cpp_torsions <corpus.smi>\n";
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "cannot open " << argv[1] << "\n";
    return 2;
  }

  std::vector<std::string> smis;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string smi;
    if (ss >> smi) {
      smis.push_back(smi);
    }
  }

  auto run_one = [](const std::string &smi, long &n_tors) -> bool {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smi));
    if (!mol) return false;
    std::unique_ptr<RDKit::ROMol> molh(RDKit::MolOps::addHs(*mol));
    ForceFields::CrystalFF::CrystalFFDetails details;
    ForceFields::CrystalFF::getExperimentalTorsions(
        *molh, details,
        /*useExpTorsions=*/true, /*useSmallRingTorsions=*/false,
        /*useMacrocycleTorsions=*/true, /*useBasicKnowledge=*/true,
        /*version=*/2, /*verbose=*/false);
    n_tors += static_cast<long>(details.expTorsionAtoms.size());
    return true;
  };

  // warm-up: loads the SMARTS torsion library and its caches
  long warm = 0;
  for (size_t i = 0; i < smis.size() && i < 5; ++i) run_one(smis[i], warm);

  long n_ok = 0, n_tors = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (const auto &smi : smis) {
    if (run_one(smi, n_tors)) ++n_ok;
  }
  auto t1 = std::chrono::steady_clock::now();
  double secs = std::chrono::duration<double>(t1 - t0).count();

  std::cerr << "cpp rdkit e2e (SMILES->AddHs->getExperimentalTorsions): " << n_ok
            << " mols in " << secs << " s = " << (1000.0 * secs / n_ok)
            << " ms/mol, " << n_tors << " torsions\n";
  return 0;
}
