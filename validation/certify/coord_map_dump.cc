// coord_map / bounds dumper — mirrors bb-rdkit/src/bridge.cc `coord_map_bounds` + `compute_raw_bounds`,
// linked against HER container's /soft/rdkit_libs. Certifies the coord-map-adjusted, triangle-smoothed
// bounds matrix (setupInitialBoundsMatrix success path). With k=0 pins it is the plain smoothed bounds
// matrix, so this one dumper certifies both "bounds" and "coord_map".
//
//   coord_map_dump <smiles>   (stdin: k  then k lines of "idx x y z")
// prints the n*n bounds matrix row-major (one f64 per line). getUpperBound above the diagonal,
// getLowerBound below — matching the bridge / native spec.bounds layout.
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <Geometry/point.h>
#include <cstdio>
#include <map>
#include <memory>

// adjustBoundsMatFromCoordMap is exported by libRDKitDistGeomHelpers with no public header —
// forward-declare exactly as bridge.cc does.
namespace RDKit {
namespace DGeomHelpers {
namespace EmbeddingOps {
void adjustBoundsMatFromCoordMap(DistGeom::BoundsMatPtr mmat, unsigned int nAtoms,
                                 const std::map<int, RDGeom::Point3D> *coordMap);
}  // namespace EmbeddingOps
}  // namespace DGeomHelpers
}  // namespace RDKit

static DistGeom::BoundsMatPtr compute_raw_bounds(const RDKit::ROMol &mol) {
  unsigned int n = mol.getNumAtoms();
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
  RDKit::DGeomHelpers::initBoundsMat(mat);
  RDKit::DGeomHelpers::setTopolBounds(mol, mat, true, false, true);
  return mat;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: coord_map_dump <smiles>   (stdin: k then k lines 'idx x y z')\n");
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

  // stdin comes from certify.py as a uniform stream of decimal doubles: k, then k*(idx x y z).
  // Read everything as double and cast the integer fields, so "3.0"/"0.0" parse cleanly.
  double kd = 0.0;
  if (scanf("%lf", &kd) != 1) {
    fprintf(stderr, "PIN_COUNT_FAIL\n");
    return 1;
  }
  int k = (int)kd;
  std::map<int, RDGeom::Point3D> coordMap;
  for (int p = 0; p < k; ++p) {
    double idxd, x, y, z;
    if (scanf("%lf %lf %lf %lf", &idxd, &x, &y, &z) != 4) {
      fprintf(stderr, "PIN_FAIL at %d\n", p);
      return 1;
    }
    coordMap[(int)idxd] = RDGeom::Point3D(x, y, z);
  }

  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  RDKit::DGeomHelpers::EmbeddingOps::adjustBoundsMatFromCoordMap(mat, n, &coordMap);
  if (!DistGeom::triangleSmoothBounds(mat, 0.05)) {
    fprintf(stderr, "SMOOTH_FAIL\n");
    return 1;
  }
  const double *data = mat->getData();
  for (std::size_t i = 0; i < (std::size_t)n * n; ++i) {
    printf("%.17g\n", data[i]);
  }
  return 0;
}
