// Stage-A/B force-field gradient dumper — mirrors bb-rdkit/src/bridge.cc `stage_a_ff_grad` +
// `compute_raw_bounds`, but linked against HER container's /soft/rdkit_libs. `constructForceField` is
// the same FF for both stages, differing only in weights: Stage-A = firstMinimization (chiral 1.0,
// fourth 0.1), Stage-B = minimizeFourthDimension (0.2, 1.0). Weights/basin are optional argv so one
// binary certifies both. Reads a SMILES (argv[1]) and n*4 coords on stdin; prints the gradient (one
// f64 per line). Diffing this against our host-bridge golden tells us whether our from-source build ≡
// her actual RDKit build.
//   stage_a_dump <smiles> [weight_chiral=1.0] [weight_fourth=0.1] [basin=5.0]   (coords n*4 on stdin)
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/Conformer.h>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/ChiralSet.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <ForceField/ForceField.h>
#include <Geometry/point.h>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>

// findChiralSets is exported by libRDKitDistGeomHelpers with no public header — forward-declare it
// exactly as bridge.cc does, so we call the real routine.
namespace RDKit {
namespace DGeomHelpers {
namespace EmbeddingOps {
void findChiralSets(const RDKit::ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters,
                    DistGeom::VECT_CHIRALSET &tetrahedralCenters,
                    const std::map<int, RDGeom::Point3D> *coordMap);
}  // namespace EmbeddingOps
}  // namespace DGeomHelpers
}  // namespace RDKit

static DistGeom::BoundsMatPtr compute_raw_bounds(const RDKit::ROMol &mol) {
  unsigned int n = mol.getNumAtoms();
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
  RDKit::DGeomHelpers::initBoundsMat(mat);
  RDKit::DGeomHelpers::setTopolBounds(mol, mat, /*set15bounds=*/true, /*scaleVDW=*/false,
                                      /*useMacrocycle14config=*/true);
  return mat;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: stage_a_dump <smiles>   (n*4 coords on stdin)\n");
    return 2;
  }
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(argv[1]));
  if (!mol) {
    fprintf(stderr, "PARSE_FAIL\n");
    return 1;
  }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  unsigned int n = mol->getNumAtoms();
  std::vector<double> coords((size_t)n * 4);
  for (size_t i = 0; i < coords.size(); ++i) {
    if (scanf("%lf", &coords[i]) != 1) {
      fprintf(stderr, "COORD_FAIL at %zu (need %u*4)\n", i, n);
      return 1;
    }
  }

  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  DistGeom::triangleSmoothBounds(mat);

  DistGeom::VECT_CHIRALSET chiralCenters, tetrahedralCenters;
  RDKit::DGeomHelpers::EmbeddingOps::findChiralSets(*mol, chiralCenters, tetrahedralCenters, nullptr);

  RDGeom::PointPtrVect positions;
  positions.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    auto *pt = new RDGeom::PointND(4);
    for (int c = 0; c < 4; ++c) {
      (*pt)[c] = coords[(size_t)i * 4 + c];
    }
    positions.push_back(pt);
  }

  // Stage-A firstMinimization weights (1.0, 0.1) by default; Stage-B minimizeFourthDimension passes
  // (0.2, 1.0). basinThresh 5.0 (BASIN_DEFAULT). All overridable via argv.
  double wChiral = argc > 2 ? atof(argv[2]) : 1.0;
  double wFourth = argc > 3 ? atof(argv[3]) : 0.1;
  double basin = argc > 4 ? atof(argv[4]) : 5.0;
  std::unique_ptr<ForceFields::ForceField> field(DistGeom::constructForceField(
      *mat, positions, chiralCenters, wChiral, wFourth, nullptr, basin, nullptr));
  field->initialize();
  std::vector<double> g((size_t)field->dimension() * n, 0.0);
  field->calcGrad(g.data());
  for (double v : g) {
    printf("%.17g\n", v);
  }
  for (auto *pt : positions) {
    delete pt;
  }
  return 0;
}
