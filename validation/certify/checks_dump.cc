// Post-embed checks dumper — mirrors bb-rdkit/src/bridge.cc `embed_checks` + `compute_raw_bounds`,
// linked against HER container's /soft/rdkit_libs. Runs RDKit's real EmbeddingOps check routines plus
// the planarity test (improper-only FF energy vs improperAtoms*0.7) on a supplied 3D conformer, and
// prints the 6 bits in native passes_checks order: tetrahedral, chiral, planarity, double-bond-
// geometry, final-chiral, double-bond-stereo (1=pass, 0=reject). Reads SMILES (argv[1]) + n*3 coords
// on stdin. The EmbedArgs struct and the EmbeddingOps functions have no public header, so they are
// forward-declared byte-for-byte from bridge.cc.
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
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <ForceField/ForceField.h>
#include <Geometry/point.h>
#include <RDGeneral/types.h>
#include <boost/dynamic_bitset.hpp>
#include <cstdio>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

// --- forward declarations (no public header), byte-for-byte from bridge.cc ---
namespace RDKit {
namespace DGeomHelpers {
namespace detail {
struct EmbedArgs {
  boost::dynamic_bitset<> *confsOk;
  bool fourD;
  INT_VECT *fragMapping;
  std::vector<std::unique_ptr<Conformer>> *confs;
  unsigned int fragIdx;
  DistGeom::BoundsMatPtr mmat;
  DistGeom::VECT_CHIRALSET const *chiralCenters;
  DistGeom::VECT_CHIRALSET const *tetrahedralCarbons;
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> const *doubleBondEnds;
  std::vector<std::pair<std::vector<unsigned int>, int>> const *stereoDoubleBonds;
  ForceFields::CrystalFF::CrystalFFDetails *etkdgDetails;
};
}  // namespace detail
namespace EmbeddingOps {
void findChiralSets(const RDKit::ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters,
                    DistGeom::VECT_CHIRALSET &tetrahedralCenters,
                    const std::map<int, RDGeom::Point3D> *coordMap);
void findDoubleBonds(
    const RDKit::ROMol &mol,
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &doubleBondEnds,
    std::vector<std::pair<std::vector<unsigned int>, int>> &stereoDoubleBonds,
    const std::map<int, RDGeom::Point3D> *coordMap);
bool checkTetrahedralCenters(const RDGeom::PointPtrVect *positions, const detail::EmbedArgs &eargs,
                             const EmbedParameters &);
bool checkChiralCenters(const RDGeom::PointPtrVect *positions, const detail::EmbedArgs &eargs,
                        const EmbedParameters &);
bool doubleBondGeometryChecks(const RDGeom::PointPtrVect &positions, const detail::EmbedArgs &eargs,
                              EmbedParameters &, double linearTol);
bool doubleBondStereoChecks(const RDGeom::PointPtrVect &positions, const detail::EmbedArgs &eargs,
                            EmbedParameters &);
bool finalChiralChecks(RDGeom::PointPtrVect *positions, const detail::EmbedArgs &eargs,
                       EmbedParameters &);
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
    fprintf(stderr, "usage: checks_dump <smiles>   (n*3 coords on stdin)\n");
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

  DistGeom::BoundsMatPtr mmat = compute_raw_bounds(*mol);
  DistGeom::triangleSmoothBounds(mmat);

  DistGeom::VECT_CHIRALSET chiralCenters, tetrahedralCarbons;
  RDKit::DGeomHelpers::EmbeddingOps::findChiralSets(*mol, chiralCenters, tetrahedralCarbons, nullptr);
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> doubleBondEnds;
  std::vector<std::pair<std::vector<unsigned int>, int>> stereoDoubleBonds;
  RDKit::DGeomHelpers::EmbeddingOps::findDoubleBonds(*mol, doubleBondEnds, stereoDoubleBonds, nullptr);

  RDGeom::PointPtrVect positions;
  positions.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    positions.push_back(new RDGeom::Point3D(coords[(size_t)i * 3 + 0], coords[(size_t)i * 3 + 1],
                                            coords[(size_t)i * 3 + 2]));
  }

  RDKit::DGeomHelpers::detail::EmbedArgs eargs;
  eargs.confsOk = nullptr;
  eargs.fourD = false;
  eargs.fragMapping = nullptr;
  eargs.confs = nullptr;
  eargs.fragIdx = 0;
  eargs.mmat = mmat;
  eargs.chiralCenters = &chiralCenters;
  eargs.tetrahedralCarbons = &tetrahedralCarbons;
  eargs.doubleBondEnds = &doubleBondEnds;
  eargs.stereoDoubleBonds = &stereoDoubleBonds;
  eargs.etkdgDetails = nullptr;

  // planarity: improper-only FF energy vs improperAtoms.size()*0.7 (transcribed from
  // minimizeWithExpTorsions, as in bridge.cc).
  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false, true, true, 2, false);
  RDGeom::Point3DPtrVect pos3d;
  pos3d.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    pos3d.push_back(new RDGeom::Point3D(coords[(size_t)i * 3 + 0], coords[(size_t)i * 3 + 1],
                                        coords[(size_t)i * 3 + 2]));
  }
  bool planar = true;
  {
    std::unique_ptr<ForceFields::ForceField> impff(
        DistGeom::construct3DImproperForceField(*mmat, pos3d, details));
    impff->initialize();
    if (impff->calcEnergy() > details.improperAtoms.size() * 0.7) {
      planar = false;
    }
  }
  for (auto *p : pos3d) {
    delete p;
  }

  RDKit::DGeomHelpers::EmbedParameters params;
  namespace EO = RDKit::DGeomHelpers::EmbeddingOps;
  int bits[6];
  bits[0] = EO::checkTetrahedralCenters(&positions, eargs, params) ? 1 : 0;
  bits[1] = EO::checkChiralCenters(&positions, eargs, params) ? 1 : 0;
  bits[2] = planar ? 1 : 0;
  bits[3] = EO::doubleBondGeometryChecks(positions, eargs, params, 1e-3) ? 1 : 0;
  bits[4] = EO::finalChiralChecks(&positions, eargs, params) ? 1 : 0;
  bits[5] = EO::doubleBondStereoChecks(positions, eargs, params) ? 1 : 0;
  for (int b = 0; b < 6; ++b) {
    printf("%d\n", bits[b]);
  }
  for (auto *p : positions) {
    delete p;
  }
  return 0;
}
