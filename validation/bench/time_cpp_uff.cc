// C++ baseline for UFF atom typing, against the same patched RDKit the bridge links.
//
//   time_cpp_uff <corpus.smi>
//
// The speed floor for bb_perceive::uff. It times the typing component in isolation: molecules are
// prepared (SmilesToMol + addHs + findSSSR) up front, outside the timer; only UFF::getAtomTypes is
// timed, after a warm-up. getAtomTypes resolves each atom to its AtomicParams (getAtomLabel + table
// lookup) — the same work the Rust bench does (atom_label + params_for_label).
//
// Equal-work proof (must match the Rust bench): atoms_typed (integer) and a checksum = sum of r1
// over all typed atoms (float). It calls RDKit's own getAtomTypes; a stopwatch, not a reimpl.
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <ForceField/UFF/Params.h>

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: time_cpp_uff <corpus.smi>\n";
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "cannot open " << argv[1] << "\n";
    return 2;
  }

  std::vector<std::unique_ptr<RDKit::ROMol>> mols;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string smi;
    if (!(ss >> smi)) continue;
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smi));
    if (!mol) continue;
    std::unique_ptr<RDKit::ROMol> molh(RDKit::MolOps::addHs(*mol));
    if (!molh->getRingInfo()->isInitialized()) {
      RDKit::MolOps::findSSSR(*molh);
    }
    mols.push_back(std::move(molh));
  }

  auto type_one = [](const RDKit::ROMol &mol, long &atoms, double &checksum) {
    auto [params, foundAll] = RDKit::UFF::getAtomTypes(mol);
    (void)foundAll;
    for (const auto *p : params) {
      if (p) {
        checksum += p->r1;
        ++atoms;
      }
    }
  };

  {
    long a = 0;
    double s = 0.0;
    for (size_t i = 0; i < mols.size() && i < 5; ++i) type_one(*mols[i], a, s);
  }

  long atoms = 0;
  double checksum = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (const auto &mol : mols) {
    type_one(*mol, atoms, checksum);
  }
  auto t1 = std::chrono::steady_clock::now();
  double secs = std::chrono::duration<double>(t1 - t0).count();

  std::cerr << "cpp rdkit UFF getAtomTypes (prepared mols): " << mols.size()
            << " mols in " << secs << " s = " << (1000.0 * secs / mols.size())
            << " ms/mol, atoms_typed " << atoms << ", checksum "
            << std::setprecision(12) << checksum << "\n";
  return 0;
}
