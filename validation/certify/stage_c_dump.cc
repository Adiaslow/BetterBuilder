// Stage-C force-field gradient dumper — mirrors bb-rdkit/src/bridge.cc `stage_c_ff_grad`, linked
// against HER container's /soft/rdkit_libs. The amide-patch-critical call is getExperimentalTorsions,
// which reads the patched torsion-preference table (her 80.0/130.9 amide force constants). Reads a
// SMILES (argv[1]) and n*3 coords on stdin; prints the gradient (one f64 per line). Diffing this vs
// our host-bridge Stage-C output certifies her amide patch is byte-faithfully reproduced by our build.
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <ForceField/ForceField.h>
#include <Geometry/point.h>
#include <cstdio>
#include <memory>
#include <vector>

static DistGeom::BoundsMatPtr compute_raw_bounds(const RDKit::ROMol &mol) {
  unsigned int n = mol.getNumAtoms();
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
  RDKit::DGeomHelpers::initBoundsMat(mat);
  RDKit::DGeomHelpers::setTopolBounds(mol, mat, true, false, true);
  return mat;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: stage_c_dump <smiles>   (n*3 coords on stdin)\n");
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
  std::vector<double> coords((size_t)n * 3);
  for (size_t i = 0; i < coords.size(); ++i) {
    if (scanf("%lf", &coords[i]) != 1) {
      fprintf(stderr, "COORD_FAIL at %zu (need %u*3)\n", i, n);
      return 1;
    }
  }

  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  DistGeom::triangleSmoothBounds(mat);

  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(
      *mol, details, /*useExpTorsions=*/true, /*useSmallRingTorsions=*/false,
      /*useMacrocycleTorsions=*/true, /*useBasicKnowledge=*/true, /*version=*/2, /*verbose=*/false);
  RDKit::DGeomHelpers::collectBondsAndAngles(*mol, details.bonds, details.angles);
  details.boundsMatForceScaling = 1.0;

  RDGeom::Point3DPtrVect positions3D;
  positions3D.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    positions3D.push_back(new RDGeom::Point3D(coords[(size_t)i * 3 + 0], coords[(size_t)i * 3 + 1],
                                              coords[(size_t)i * 3 + 2]));
  }

  std::unique_ptr<ForceFields::ForceField> field(
      DistGeom::construct3DForceField(*mat, positions3D, details));
  field->initialize();
  std::vector<double> g((size_t)3 * n, 0.0);
  field->calcGrad(g.data());
  for (double v : g) {
    printf("%.17g\n", v);
  }
  for (auto *pt : positions3D) {
    delete pt;
  }
  return 0;
}
