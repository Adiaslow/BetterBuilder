// C++ baseline for the bounds-matrix construction, against the same patched RDKit the bridge links.
//
//   time_cpp_bounds <corpus.smi>
//
// The speed floor for the pure-Rust setTopolBounds port. It times the component in isolation:
// molecules are read and prepared (SmilesToMol + addHs + findSSSR) up front, outside the timer;
// only initBoundsMat + setTopolBounds is timed, after a warm-up, so process startup, I/O, and
// SMILES perception are all excluded. The Rust counterpart must time the same scope.
//
// It prints two figures that MUST match the Rust bench, or the two are not doing equal work:
//   * constrained_pairs — integer count of upper bounds tightened below MAX_UPPER (exact, no float)
//   * checksum          — sum of every matrix entry (finer, float)
// It calls RDKit's own setTopolBounds; it is a stopwatch around the real library, not a reimpl.
#include <chrono>
#include <cmath>
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
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <DistGeom/BoundsMatrix.h>

// MAX_UPPER from BoundsMatrixBuilder.cpp: the default upper bound, i.e. "unconstrained".
static constexpr double MAX_UPPER = 1000.0;

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: time_cpp_bounds <corpus.smi>\n";
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "cannot open " << argv[1] << "\n";
    return 2;
  }

  // Prepare every molecule up front (outside the timer): SMILES -> AddHs -> ring info. This mirrors
  // build_mol_spec / raw_bounds, so setTopolBounds sees exactly the molecule the oracle fixture was
  // captured from.
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

  auto build_one = [](const RDKit::ROMol &mol, long &constrained, double &checksum) {
    unsigned int n = mol.getNumAtoms();
    DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
    RDKit::DGeomHelpers::initBoundsMat(mat);
    RDKit::DGeomHelpers::setTopolBounds(mol, mat, /*set15bounds=*/true,
                                        /*scaleVDW=*/false,
                                        /*useMacrocycle14config=*/true);
    const double *d = mat->getData();
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = i + 1; j < n; ++j) {
        if (d[i * n + j] < MAX_UPPER) ++constrained;  // upper bound was tightened
      }
    }
    for (std::size_t i = 0; i < static_cast<std::size_t>(n) * n; ++i) {
      checksum += d[i];
    }
  };

  // warm-up
  {
    long c = 0;
    double s = 0.0;
    for (size_t i = 0; i < mols.size() && i < 5; ++i) build_one(*mols[i], c, s);
  }

  long constrained = 0;
  double checksum = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (const auto &mol : mols) {
    build_one(*mol, constrained, checksum);
  }
  auto t1 = std::chrono::steady_clock::now();
  double secs = std::chrono::duration<double>(t1 - t0).count();

  std::cerr << "cpp rdkit setTopolBounds (prepared mols): " << mols.size()
            << " mols in " << secs << " s = " << (1000.0 * secs / mols.size())
            << " ms/mol, constrained_pairs " << constrained << ", checksum "
            << std::setprecision(12) << checksum << "\n";
  return 0;
}
