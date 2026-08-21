#include "bb-rdkit/src/bridge.h"

#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Geometry/point.h>

#include <boost/tuple/tuple.hpp>

#include <ForceField/UFF/Utils.h>
#include <ForceField/UFF/BondStretch.h>
#include <ForceField/UFF/Params.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/PeriodicTable.h>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/ChiralSet.h>
#include <DistGeom/TriangleSmooth.h>
#include <DistGeom/DistGeomUtils.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleContribs.h>
#include <GraphMol/Atom.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <RDGeneral/types.h>
#include <boost/dynamic_bitset.hpp>

// Forward-declare the oracle's exported findDoubleBonds (in EmbeddingOps; no public header) so we
// call the real routine rather than reimplementing it. Symbol lives in libRDKitDistGeomHelpers.
namespace RDKit {
namespace DGeomHelpers {
// detail::EmbedArgs is defined only in Embedder.cpp (no public header); the exported check functions
// take it by reference. Redeclared here byte-for-byte from that definition so the layout matches and
// we can hand the checks the exact struct they expect (only mmat + the chiral/double-bond sets are
// read by the checks we call; the rest stay null/default).
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
void adjustBoundsMatFromCoordMap(DistGeom::BoundsMatPtr mmat, unsigned int nAtoms,
                                 const std::map<int, RDGeom::Point3D> *coordMap);
bool checkTetrahedralCenters(const RDGeom::PointPtrVect *positions,
                             const detail::EmbedArgs &eargs, const EmbedParameters &);
bool checkChiralCenters(const RDGeom::PointPtrVect *positions,
                        const detail::EmbedArgs &eargs, const EmbedParameters &);
bool doubleBondGeometryChecks(const RDGeom::PointPtrVect &positions,
                              const detail::EmbedArgs &eargs, EmbedParameters &,
                              double linearTol);
bool doubleBondStereoChecks(const RDGeom::PointPtrVect &positions,
                            const detail::EmbedArgs &eargs, EmbedParameters &);
bool finalChiralChecks(RDGeom::PointPtrVect *positions, const detail::EmbedArgs &eargs,
                       EmbedParameters &);
}  // namespace EmbeddingOps
}  // namespace DGeomHelpers
}  // namespace RDKit

// The UFF parameter table is an `extern const std::string` in Params.cpp (external linkage via its
// extern declaration there). Forward-declare it at global scope so uff_param_table can hand the
// exact bytes to Rust rather than transcribing the table.
namespace ForceFields {
namespace UFF {
extern const std::string defaultParamData;
}  // namespace UFF
}  // namespace ForceFields

namespace bb {

std::size_t num_heavy_atoms(rust::Str smiles) {
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  return mol ? static_cast<std::size_t>(mol->getNumHeavyAtoms()) : 0;
}

// The raw (pre-triangleSmooth) bounds matrix, with exactly the flags the Python
// GetMoleculeBoundsMatrix uses (set15bounds, useMacrocycle14config; the remaining flags take their
// header defaults forceTransAmides=true, set14bounds=true, set13bounds=true). Single-sourced here so
// build_mol_spec and the raw_bounds oracle can never disagree on how the matrix is built.
static DistGeom::BoundsMatPtr compute_raw_bounds(const RDKit::ROMol &mol) {
  unsigned int n = mol.getNumAtoms();
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
  RDKit::DGeomHelpers::initBoundsMat(mat);
  RDKit::DGeomHelpers::setTopolBounds(mol, mat, /*set15bounds=*/true,
                                      /*scaleVDW=*/false,
                                      /*useMacrocycle14config=*/true);
  return mat;
}


// RDKit's Stage-A distance-geometry force field (firstMinimization's
// constructForceField(mmat, pos, csets, 1.0, 0.1, nullptr, 5.0)) evaluated at the given 4D
// coordinates. The oracle for bb-embed's stage_a_energy: same smoothed bounds + chiral sets that
// build_mol_spec emits, so any energy difference at identical coords is a force-field *formula*
// difference, not an input difference. `coords` is n*4 doubles (x,y,z,w per atom, addHs order).
// NaN signals a parse failure or a coords-length mismatch.

// RDKit's Stage-C 3D force field (construct3DForceField: experimental torsions + UFF impropers +
// 1-3 and long-range distance constraints) gradient at the given 3D coordinates. Same smoothed
// bounds and getExperimentalTorsions details build_mol_spec emits, so a gradient gap is a
// force-field formula difference. `coords` is n*3 doubles (addHs order). Empty on failure.

// RDKit's improper-only Stage-C force field (construct3DImproperForceField — the planarity-check FF)
// gradient at the given 3D coordinates. Isolates the UFF inversion term from torsions/constraints.

// RDKit's experimental-torsion term alone (TorsionAngleContribs on a bare force field) gradient at
// the given 3D coordinates — isolates the M6 torsion physics from the impropers and constraints.
rust::Vec<double> stage_c_torsion_grad(rust::Str smiles, rust::Slice<const double> coords) {
  rust::Vec<double> empty;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) { return empty; }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) { RDKit::MolOps::findSSSR(*mol); }
  unsigned int n = mol->getNumAtoms();
  if (coords.size() != static_cast<std::size_t>(n) * 3) { return empty; }
  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false, true, true, 2, false);

  ForceFields::ForceField ff(3);
  std::vector<RDGeom::Point3D *> pts;
  pts.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    auto *p = new RDGeom::Point3D(coords[(std::size_t)i*3+0], coords[(std::size_t)i*3+1], coords[(std::size_t)i*3+2]);
    pts.push_back(p);
    ff.positions().push_back(p);
  }
  auto tc = std::make_unique<ForceFields::CrystalFF::TorsionAngleContribs>(&ff);
  for (std::size_t t = 0; t < details.expTorsionAtoms.size(); ++t) {
    tc->addContrib(details.expTorsionAtoms[t][0], details.expTorsionAtoms[t][1],
                   details.expTorsionAtoms[t][2], details.expTorsionAtoms[t][3],
                   details.expTorsionAngles[t].second, details.expTorsionAngles[t].first);
  }
  rust::Vec<double> grad;
  if (!tc->empty()) {
    ff.contribs().push_back(std::move(tc));
  }
  ff.initialize();
  std::vector<double> g((std::size_t)3 * n, 0.0);
  ff.calcGrad(g.data());
  grad.reserve(g.size());
  for (double v : g) grad.push_back(v);
  for (auto *p : pts) delete p;
  return grad;
}

rust::Vec<double> stage_c_improper_grad(rust::Str smiles, rust::Slice<const double> coords) {
  rust::Vec<double> empty;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) { return empty; }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) { RDKit::MolOps::findSSSR(*mol); }
  unsigned int n = mol->getNumAtoms();
  if (coords.size() != static_cast<std::size_t>(n) * 3) { return empty; }
  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  DistGeom::triangleSmoothBounds(mat);
  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(
      *mol, details, true, false, true, true, 2, false);
  RDGeom::Point3DPtrVect positions3D;
  positions3D.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    positions3D.push_back(new RDGeom::Point3D(coords[(std::size_t)i*3+0], coords[(std::size_t)i*3+1], coords[(std::size_t)i*3+2]));
  }
  rust::Vec<double> grad;
  {
    std::unique_ptr<ForceFields::ForceField> field(
        DistGeom::construct3DImproperForceField(*mat, positions3D, details));
    field->initialize();
    std::vector<double> g((std::size_t)3 * n, 0.0);
    field->calcGrad(g.data());
    grad.reserve(g.size());
    for (double v : g) grad.push_back(v);
  }
  for (auto *pt : positions3D) delete pt;
  return grad;
}

// RDKit's post-embed accept/reject checks (EmbeddingOps::check*) run on the given 3D conformer, in
// the order the embed loop applies them. Returns 5 bits {tetrahedral, chiral, double-bond-geometry,
// final-chiral, double-bond-stereo}. The chiral/tetrahedral sets and double-bond ends come from the
// same exported findChiralSets/findDoubleBonds the embedder uses; mmat is the smoothed bounds matrix
// (finalChiralChecks' distance test reads it). This is the oracle for bb_embed::checks.
rust::Vec<int> embed_checks(rust::Str smiles, rust::Slice<const double> coords) {
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) { return out; }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) { RDKit::MolOps::findSSSR(*mol); }
  unsigned int n = mol->getNumAtoms();
  if (coords.size() != static_cast<std::size_t>(n) * 3) { return out; }

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
    positions.push_back(new RDGeom::Point3D(coords[(std::size_t)i * 3 + 0],
                                            coords[(std::size_t)i * 3 + 1],
                                            coords[(std::size_t)i * 3 + 2]));
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

  // The planarity check isn't a standalone exported function — it lives inside minimizeWithExpTorsions:
  // an improper-only FF whose energy must not exceed improperAtoms.size()*0.7. Transcribed here (the
  // etkdgDetails.improperAtoms come from getExperimentalTorsions, as in the embedder).
  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false, true, true, 2, false);
  RDGeom::Point3DPtrVect pos3d;
  pos3d.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    pos3d.push_back(new RDGeom::Point3D(coords[(std::size_t)i * 3 + 0],
                                        coords[(std::size_t)i * 3 + 1],
                                        coords[(std::size_t)i * 3 + 2]));
  }
  bool planar = true;
  {
    std::unique_ptr<ForceFields::ForceField> impff(
        DistGeom::construct3DImproperForceField(*mmat, pos3d, details));
    impff->initialize();
    if (impff->calcEnergy() > details.improperAtoms.size() * 0.7) { planar = false; }
  }
  for (auto *p : pos3d) delete p;

  RDKit::DGeomHelpers::EmbedParameters params;  // defaults; trackFailures=false
  namespace EO = RDKit::DGeomHelpers::EmbeddingOps;
  // Order matches native passes_checks: tetrahedral, chiral, planarity, double-bond-geometry,
  // final-chiral, double-bond-stereo.
  out.push_back(EO::checkTetrahedralCenters(&positions, eargs, params) ? 1 : 0);
  out.push_back(EO::checkChiralCenters(&positions, eargs, params) ? 1 : 0);
  out.push_back(planar ? 1 : 0);
  out.push_back(EO::doubleBondGeometryChecks(positions, eargs, params, 1e-3) ? 1 : 0);
  out.push_back(EO::finalChiralChecks(&positions, eargs, params) ? 1 : 0);
  out.push_back(EO::doubleBondStereoChecks(positions, eargs, params) ? 1 : 0);

  for (auto *p : positions) delete p;
  return out;
}

rust::Vec<double> stage_c_ff_grad(rust::Str smiles, rust::Slice<const double> coords) {
  rust::Vec<double> empty;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return empty;
  }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  unsigned int n = mol->getNumAtoms();
  if (coords.size() != static_cast<std::size_t>(n) * 3) {
    return empty;
  }
  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  DistGeom::triangleSmoothBounds(mat);

  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(
      *mol, details, /*useExpTorsions=*/true, /*useSmallRingTorsions=*/false,
      /*useMacrocycleTorsions=*/true, /*useBasicKnowledge=*/true, /*version=*/2, /*verbose=*/false);
  // getExperimentalTorsions fills only the torsion/improper details; the 1-2 bonds and 1-3 angles
  // that add12Terms/add13Terms pin come from collectBondsAndAngles (the embedder does this before
  // construct3DForceField). Without it every pair falls through to addLongRangeDistanceConstraints.
  RDKit::DGeomHelpers::collectBondsAndAngles(*mol, details.bonds, details.angles);
  // The embedder copies params.boundsMatForceScaling into etkdgDetails (Embedder.cpp:1252);
  // ETKDGv3 sets it to 1. addLongRangeDistanceConstraints scales its force constant by this
  // (fc = boundsMatForceScaling*10), so leaving it 0 zeroes every long-range term.
  details.boundsMatForceScaling = 1.0;

  RDGeom::Point3DPtrVect positions3D;
  positions3D.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    positions3D.push_back(new RDGeom::Point3D(coords[static_cast<std::size_t>(i) * 3 + 0],
                                              coords[static_cast<std::size_t>(i) * 3 + 1],
                                              coords[static_cast<std::size_t>(i) * 3 + 2]));
  }

  rust::Vec<double> grad;
  {
    std::unique_ptr<ForceFields::ForceField> field(
        DistGeom::construct3DForceField(*mat, positions3D, details));
    field->initialize();
    std::vector<double> g(static_cast<std::size_t>(3) * n, 0.0);
    field->calcGrad(g.data());
    grad.reserve(g.size());
    for (double v : g) {
      grad.push_back(v);
    }
  }
  for (auto *pt : positions3D) {
    delete pt;
  }
  return grad;
}

rust::Vec<double> stage_a_ff_grad(rust::Str smiles, rust::Slice<const double> coords,
                                  double weight_chiral, double weight_fourth, double basin) {
  rust::Vec<double> empty;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return empty;
  }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  unsigned int n = mol->getNumAtoms();
  if (coords.size() != static_cast<std::size_t>(n) * 4) {
    return empty;
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
      (*pt)[c] = coords[static_cast<std::size_t>(i) * 4 + c];
    }
    positions.push_back(pt);
  }

  rust::Vec<double> grad;
  {
    std::unique_ptr<ForceFields::ForceField> field(DistGeom::constructForceField(
        *mat, positions, chiralCenters, weight_chiral, weight_fourth, nullptr, basin, nullptr));
    field->initialize();
    // The force field is what drives the embed; compare gradients (forces), not the energy value —
    // RDKit's Chiral/FourthDim contribs use an energy that is inconsistent with their own gradient
    // (energy 2x the gradient's integral), so only the gradient is a meaningful cross-check.
    std::vector<double> g(static_cast<std::size_t>(field->dimension()) * n, 0.0);
    field->calcGrad(g.data());
    grad.reserve(g.size());
    for (double v : g) {
      grad.push_back(v);
    }
  }
  for (auto *pt : positions) {
    delete pt;
  }
  return grad;
}

std::unique_ptr<MolSpecC> build_mol_spec(rust::Str smiles) {
  auto out = std::make_unique<MolSpecC>();
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;
  }

  // --- two-stage recipe conformer counts: count_exo_rotatable on the HEAVY-atom graph (pre-addHs),
  //     matching build_ligands.py (exocyclic single bonds touching a ring). ---
  {
    if (!mol->getRingInfo()->isInitialized()) {
      RDKit::MolOps::findSSSR(*mol);
    }
    const RDKit::RingInfo *ri0 = mol->getRingInfo();
    std::set<int> ring_atoms;
    for (const auto &ring : ri0->atomRings()) {
      for (int a : ring) {
        ring_atoms.insert(a);
      }
    }
    int exo = 0;
    for (const auto *bond : mol->bonds()) {
      if (ri0->numBondRings(bond->getIdx()) > 0) {
        continue; // skip ring bonds
      }
      if (bond->getBondTypeAsDouble() == 1.0) {
        int i = bond->getBeginAtomIdx();
        int j = bond->getEndAtomIdx();
        if (ring_atoms.count(i) || ring_atoms.count(j)) {
          ++exo;
        }
      }
    }
    out->core_seeds_ = (exo <= 2) ? 20 : 10;
    out->sidechain_confs_ = (exo <= 2) ? 10 : 20;
  }

  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  unsigned int n = mol->getNumAtoms();
  out->n_atoms_ = n;

  // --- per-atom element + net formal charge. Taken after addHs so the indices match every other
  //     array here, and the order AMSOL reports its per-atom results in. ---
  out->atomic_nums_.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    out->atomic_nums_.push_back(
        static_cast<std::uint8_t>(mol->getAtomWithIdx(i)->getAtomicNum()));
  }
  out->formal_charge_ = static_cast<std::int32_t>(RDKit::MolOps::getFormalCharge(*mol));

  // --- bounds matrix: exactly as the Python GetMoleculeBoundsMatrix (doTriangleSmoothing,
  //     useMacrocycle14config) ---
  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  DistGeom::triangleSmoothBounds(mat);
  const double *data = mat->getData();
  out->bounds_.assign(data, data + static_cast<std::size_t>(n) * n);

  // --- experimental torsions (ETKDGv3 flags; carries Divya's patched amide reweighting) ---
  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(
      *mol, details, /*useExpTorsions=*/true, /*useSmallRingTorsions=*/false,
      /*useMacrocycleTorsions=*/true, /*useBasicKnowledge=*/true, /*version=*/2,
      /*verbose=*/false);
  for (std::size_t i = 0; i < details.expTorsionAtoms.size(); ++i) {
    for (int a : details.expTorsionAtoms[i]) {
      out->tors_atoms_.push_back(static_cast<std::uint32_t>(a));
    }
    for (int s : details.expTorsionAngles[i].first) {
      out->tors_signs_.push_back(s);
    }
    for (double v : details.expTorsionAngles[i].second) {
      out->tors_v_.push_back(v);
    }
  }

  // --- chiral sets: RDKit's own findChiralSets, forward-declared above (it has no public header
  //     but is an exported symbol). vol (0,0) == tetrahedral centre, used for checks only. ---
  {
    DistGeom::VECT_CHIRALSET chiralCenters, tetrahedralCenters;
    RDKit::DGeomHelpers::EmbeddingOps::findChiralSets(*mol, chiralCenters, tetrahedralCenters,
                                                      nullptr);
    auto emit = [&](const DistGeom::VECT_CHIRALSET &sets) {
      for (const auto &cs : sets) {
        out->chi_center_.push_back(static_cast<std::uint32_t>(cs->d_idx0));
        out->chi_atoms_.push_back(static_cast<std::uint32_t>(cs->d_idx1));
        out->chi_atoms_.push_back(static_cast<std::uint32_t>(cs->d_idx2));
        out->chi_atoms_.push_back(static_cast<std::uint32_t>(cs->d_idx3));
        out->chi_atoms_.push_back(static_cast<std::uint32_t>(cs->d_idx4));
        out->chi_vol_.push_back(cs->d_volumeLowerBound);
        out->chi_vol_.push_back(cs->d_volumeUpperBound);
        out->chi_fused_.push_back(static_cast<std::uint8_t>(
            (cs->d_structureFlags &
             static_cast<std::uint64_t>(DistGeom::ChiralSetStructureFlags::IN_FUSED_SMALL_RINGS))
                ? 1
                : 0));
      }
    };
    emit(chiralCenters);
    emit(tetrahedralCenters);
  }


  // --- UFF impropers: details.improperAtoms was populated by getExperimentalTorsions
  //     (useBasicKnowledge=true). Reproduce addImproperTorsionTerms: 3 permutation contribs per
  //     sp2 center, coeffs from RDKit's own calcInversionCoefficientsAndForceConstant, fc = res*10
  //     (kTermImproper). improperAtoms entry = [nbr0, center, nbr1, nbr2, atomicNum, isCBoundToO]. ---
  static const int perm[3][3] = {{0, 2, 3}, {0, 3, 2}, {2, 3, 0}}; // (n0, n2, n3); center is idx1
  for (const auto &imp : details.improperAtoms) {
    int atomicNum = imp[4];
    bool isCBoundToO = static_cast<bool>(imp[5]);
    auto [res, C0, C1, C2] =
        ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant(
            atomicNum, isCBoundToO);
    double fc = res * 10.0;
    for (int p = 0; p < 3; ++p) {
      out->imp_atoms_.push_back(static_cast<std::uint32_t>(imp[perm[p][0]]));
      out->imp_atoms_.push_back(static_cast<std::uint32_t>(imp[1])); // center
      out->imp_atoms_.push_back(static_cast<std::uint32_t>(imp[perm[p][1]]));
      out->imp_atoms_.push_back(static_cast<std::uint32_t>(imp[perm[p][2]]));
      out->imp_coef_.push_back(C0);
      out->imp_coef_.push_back(C1);
      out->imp_coef_.push_back(C2);
      out->imp_coef_.push_back(fc);
    }
  }

  // --- bonds & angles (for Stage-C 1-2 / 1-3 distance constraints). collectBondsAndAngles is the
  //     exact routine Embedder.cpp feeds into construct3DForceField. angle = [i, j=center, k, flag]
  //     where flag=1 for triple bonds / consecutive doubles (→ 179-180° angle constraint). ---
  std::vector<std::pair<int, int>> bonds;
  std::vector<std::vector<int>> angles;
  RDKit::DGeomHelpers::collectBondsAndAngles(*mol, bonds, angles);
  for (const auto &b : bonds) {
    out->bonds_.push_back(static_cast<std::uint32_t>(b.first));
    out->bonds_.push_back(static_cast<std::uint32_t>(b.second));
  }
  for (const auto &a : angles) {
    out->angles_.push_back(static_cast<std::uint32_t>(a[0]));
    out->angles_.push_back(static_cast<std::uint32_t>(a[1]));
    out->angles_.push_back(static_cast<std::uint32_t>(a[2]));
    out->angles_.push_back(static_cast<std::uint32_t>(a[3]));
  }

  // --- core-pin recipe atoms: largest symmetrized-SSSR ring (GetSymmSSSR) + exocyclic pins.
  //     Replays build_ligands.py's cmap construction: every core atom, plus exocyclic =O and
  //     exocyclic C bonded to a ring N. bb-embed freezes these to the seed conformer's coords. ---
  {
    std::vector<std::vector<int>> symmRings;
    RDKit::MolOps::symmetrizeSSSR(*mol, symmRings);
    if (!symmRings.empty()) {
      const std::vector<int> *largest = &symmRings[0];
      for (const auto &r : symmRings) {
        if (r.size() > largest->size()) {
          largest = &r;
        }
      }
      std::set<int> core(largest->begin(), largest->end());
      std::set<int> pin(core.begin(), core.end());
      for (int k : core) {
        const RDKit::Atom *ak = mol->getAtomWithIdx(k);
        RDKit::ROMol::OEDGE_ITER b, e;
        boost::tie(b, e) = mol->getAtomBonds(ak);
        while (b != e) {
          const RDKit::Bond *bond = (*mol)[*b];
          const RDKit::Atom *nbr = bond->getOtherAtom(ak);
          int nb = nbr->getIdx();
          // exocyclic double-bonded O (e.g. carbonyl)
          if (nbr->getAtomicNum() == 8 && bond->getBondTypeAsDouble() == 2.0 &&
              !pin.count(nb)) {
            pin.insert(nb);
          }
          // exocyclic C on a ring N
          if (ak->getAtomicNum() == 7 && nbr->getAtomicNum() == 6 &&
              !core.count(nb) && !pin.count(nb)) {
            pin.insert(nb);
          }
          ++b;
        }
      }
      for (int a : pin) {
        out->pin_atoms_.push_back(static_cast<std::uint32_t>(a));
      }
    }
  }

  // --- double bonds (linear + stereo acceptance checks), via the oracle's own findDoubleBonds
  //     (coordMap nullptr: extract all; pinned stereo bonds simply pass the check trivially). ---
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> dbEnds;
  std::vector<std::pair<std::vector<unsigned int>, int>> stereoDB;
  RDKit::DGeomHelpers::EmbeddingOps::findDoubleBonds(*mol, dbEnds, stereoDB, nullptr);
  for (const auto &t : dbEnds) {
    out->dbe_.push_back(static_cast<std::uint32_t>(std::get<0>(t)));
    out->dbe_.push_back(static_cast<std::uint32_t>(std::get<1>(t)));
    out->dbe_.push_back(static_cast<std::uint32_t>(std::get<2>(t)));
  }
  for (const auto &s : stereoDB) {
    for (unsigned int a : s.first) {
      out->sdb_atoms_.push_back(static_cast<std::uint32_t>(a));
    }
    out->sdb_sign_.push_back(s.second);
  }

  out->ok_ = true;
  return out;
}

rust::String uff_param_table() {
  return rust::String(::ForceFields::UFF::defaultParamData);
}

rust::Vec<double> raw_bounds_stage(rust::Str smiles, bool set13, bool set14, bool set15) {
  rust::Vec<double> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;  // empty: RDKit could not parse this SMILES
  }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  // Same molecule and same non-staged flags as compute_raw_bounds (scaleVDW=false,
  // useMacrocycle14config=true, forceTransAmides=true), but with set13/set14/set15 selectable so the
  // bounds-matrix port can be gated one stage at a time (e.g. set13=set14=set15=false is set12+VDW).
  unsigned int n = mol->getNumAtoms();
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
  RDKit::DGeomHelpers::initBoundsMat(mat);
  RDKit::DGeomHelpers::setTopolBounds(*mol, mat, /*set15bounds=*/set15,
                                      /*scaleVDW=*/false,
                                      /*useMacrocycle14config=*/true,
                                      /*forceTransAmides=*/true,
                                      /*set14bounds=*/set14, /*set13bounds=*/set13);
  const double *data = mat->getData();
  out.reserve(static_cast<std::size_t>(n) * n);
  for (std::size_t i = 0; i < static_cast<std::size_t>(n) * n; ++i) {
    out.push_back(data[i]);
  }
  return out;
}

rust::Vec<int> chiral_tags(rust::Str smiles) {
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;  // empty: RDKit could not parse this SMILES
  }
  RDKit::MolOps::addHs(*mol);
  // per-atom getChiralTag() in post-AddHs order: CHI_UNSPECIFIED=0, CHI_TETRAHEDRAL_CW=1,
  // CHI_TETRAHEDRAL_CCW=2, ... — the oracle for the Rust tetrahedral-chirality port.
  for (const auto atom : mol->atoms()) {
    out.push_back(static_cast<int>(atom->getChiralTag()));
  }
  return out;
}

rust::Vec<int> set_arom_only(rust::Str smiles) {
  // Isolate MolOps::setAromaticity: parse, updatePropertyCache, Kekulize, setAromaticity — and stop.
  // Per-atom isAromatic in SMILES order; leading -1 = failure. Tells whether a LATER sanitize step
  // (not setAromaticity itself) is responsible for a difference.
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles), 0, false));
  if (!mol) { out.push_back(-1); return out; }
  try {
    mol->updatePropertyCache();
    RDKit::MolOps::Kekulize(*mol, true);
    RDKit::MolOps::setAromaticity(*mol);
  } catch (...) { out.push_back(-1); return out; }
  for (const auto atom : mol->atoms()) out.push_back(atom->getIsAromatic() ? 1 : 0);
  return out;
}

rust::Vec<int> atom_elec(rust::Str smiles) {
  // RDKit MolOps::countAtomElec per atom, on the KEKULIZED molecule — the state setAromaticity sees.
  // Heavy-atom (SMILES) order. Leading -1 = parse failure.
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles), 0, false));
  if (!mol) {
    out.push_back(-1);
    return out;
  }
  try {
    mol->updatePropertyCache();
    RDKit::MolOps::Kekulize(*mol, true);
  } catch (...) {
    out.push_back(-1);
    return out;
  }
  for (const auto atom : mol->atoms()) {
    out.push_back(RDKit::MolOps::countAtomElec(atom));
  }
  return out;
}

rust::Vec<int> aromatic_perception(rust::Str smiles) {
  // RDKit's PERCEIVED aromaticity after sanitization, on the heavy-atom mol (SMILES atom order,
  // pre-AddHs). This is the ground truth for bb-perceive's aromaticity port: RDKit re-perceives
  // aromaticity regardless of how the input wrote it (lowercase vs Kekulé), which bb-perceive must
  // reproduce. Layout: first one int per atom (isAromatic 0/1), then 3 ints per bond
  // (beginIdx, endIdx, isAromatic). A leading -1 marks a parse failure (empty otherwise is ambiguous
  // with a 0-atom mol).
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    out.push_back(-1);
    return out;
  }
  out.push_back(static_cast<int>(mol->getNumAtoms()));
  for (const auto atom : mol->atoms()) {
    out.push_back(atom->getIsAromatic() ? 1 : 0);
  }
  for (const auto bond : mol->bonds()) {
    out.push_back(static_cast<int>(bond->getBeginAtomIdx()));
    out.push_back(static_cast<int>(bond->getEndAtomIdx()));
    out.push_back(bond->getIsAromatic() ? 1 : 0);
  }
  return out;
}

rust::Vec<int> bond_dirs(rust::Str smiles) {
  // Per bond, AS PARSED (sanitize=false, so bond directions survive): (beginIdx, endIdx, bondDir).
  // BondDir: NONE=0, BEGINWEDGE=1, BEGINDASH=2, ENDDOWNRIGHT=3, ENDUPRIGHT=4. Ground truth for how
  // RDKit orients ring-closure directional bonds. Leading -1 marks a parse failure.
  rust::Vec<int> out;
  RDKit::SmilesParserParams ps;
  ps.sanitize = false;
  ps.removeHs = false;
  std::unique_ptr<RDKit::RWMol> mol;
  try {
    mol.reset(RDKit::SmilesToMol(std::string(smiles), ps));
  } catch (...) {
    mol.reset();
  }
  if (!mol) {
    out.push_back(-1);
    return out;
  }
  for (const auto bond : mol->bonds()) {
    out.push_back(static_cast<int>(bond->getBeginAtomIdx()));
    out.push_back(static_cast<int>(bond->getEndAtomIdx()));
    out.push_back(static_cast<int>(bond->getBondDir()));
  }
  return out;
}

rust::Vec<int> bond_types(rust::Str smiles) {
  // Per bond after sanitize: (beginIdx, endIdx, bondType enum). RDKit BondType: SINGLE=1, DOUBLE=2,
  // TRIPLE=3, AROMATIC=12, DATIVE=17 (see Bond.h). Ground truth for whether e.g. ligand→metal bonds
  // became dative. Leading -1 marks a parse failure.
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    out.push_back(-1);
    return out;
  }
  for (const auto bond : mol->bonds()) {
    out.push_back(static_cast<int>(bond->getBeginAtomIdx()));
    out.push_back(static_cast<int>(bond->getEndAtomIdx()));
    out.push_back(static_cast<int>(bond->getBondType()));
  }
  return out;
}

rust::Vec<int> sssr_rings(rust::Str smiles) {
  // RDKit's symmetrized SSSR (getRingInfo()->atomRings() after sanitize) — the ground truth for
  // bb-perceive's ring port. Flat encoding: each ring's atom indices (SMILES/heavy-atom order)
  // followed by a -1 separator. A single leading -1 with nothing after marks a parse failure.
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    out.push_back(-1);
    return out;
  }
  const auto &rings = mol->getRingInfo()->atomRings();
  for (const auto &ring : rings) {
    for (int a : ring) {
      out.push_back(a);
    }
    out.push_back(-1);
  }
  return out;
}

rust::Vec<int> bond_stereo(rust::Str smiles) {
  rust::Vec<int> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;  // empty: RDKit could not parse this SMILES
  }
  RDKit::MolOps::addHs(*mol);
  // getStereo() and getStereoAtoms() as SmilesToMol assigned them (via
  // setBondStereoFromDirections). Three ints per bond in post-AddHs bond order:
  // stereo enum (NONE=0 ANY=1 Z=2 E=3 CIS=4 TRANS=5), then the two stereo reference atoms (-1 if
  // fewer than two are set). The oracle for the Rust double-bond-stereo port.
  for (const auto bond : mol->bonds()) {
    out.push_back(static_cast<int>(bond->getStereo()));
    const auto &sa = bond->getStereoAtoms();
    out.push_back(sa.size() >= 1 ? sa[0] : -1);
    out.push_back(sa.size() >= 2 ? sa[1] : -1);
  }
  return out;
}

rust::Vec<double> topo_distance_matrix(rust::Str smiles) {
  rust::Vec<double> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;  // empty: RDKit could not parse this SMILES
  }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  // getDistanceMat with setTopolBounds' defaults (useBO=false, useAtomWts=false); the returned
  // pointer is cached on the molecule, owned by it, so it is copied out, not freed here.
  std::size_t n = mol->getNumAtoms();
  const double *d = RDKit::MolOps::getDistanceMat(*mol, false, false, false);
  out.reserve(n * n);
  for (std::size_t i = 0; i < n * n; ++i) {
    out.push_back(d[i]);
  }
  return out;
}

rust::Vec<double> uff_bond_rest_lengths(rust::Str smiles) {
  rust::Vec<double> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;  // empty: RDKit could not parse this SMILES
  }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  // Exactly what set12Bounds stores in accumData.bondLengths, per bond in RDKit bond order:
  // calcBondRestLength from the two atoms' UFF params, or the crude (vdw1+vdw2)/2 when either atom
  // has no params. bond order is getBondTypeAsDouble (1.5 for aromatic).
  auto [params, foundAll] = RDKit::UFF::getAtomTypes(*mol);
  (void)foundAll;
  const auto *ptable = RDKit::PeriodicTable::getTable();
  for (const auto bond : mol->bonds()) {
    auto b = bond->getBeginAtomIdx();
    auto e = bond->getEndAtomIdx();
    double bo = bond->getBondTypeAsDouble();
    double bl;
    if (params[b] && params[e] && bo > 0) {
      bl = ForceFields::UFF::Utils::calcBondRestLength(bo, params[b], params[e]);
    } else {
      double vw1 = ptable->getRvdw(mol->getAtomWithIdx(b)->getAtomicNum());
      double vw2 = ptable->getRvdw(mol->getAtomWithIdx(e)->getAtomicNum());
      bl = (vw1 + vw2) / 2;
    }
    out.push_back(bl);
  }
  return out;
}

rust::Vec<rust::String> uff_atom_labels(rust::Str smiles) {
  rust::Vec<rust::String> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;  // empty: RDKit could not parse this SMILES
  }
  // Same preparation as raw_bounds / build_mol_spec: getAtomTypes (hence getAtomLabel) runs on the
  // post-AddHs molecule inside set12Bounds, so labels must be taken in that same state and order.
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  out.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    out.push_back(rust::String(RDKit::UFF::Tools::getAtomLabel(atom)));
  }
  return out;
}

rust::Vec<double> raw_bounds(rust::Str smiles) {
  rust::Vec<double> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    return out;  // empty: RDKit could not parse this SMILES
  }
  // Same preparation as build_mol_spec: hydrogens added, then ring info ensured. setTopolBounds
  // reads ring info throughout, so it must be initialized before the call.
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  std::size_t n = mol->getNumAtoms();
  const double *data = mat->getData();
  out.reserve(n * n);
  for (std::size_t i = 0; i < n * n; ++i) {
    out.push_back(data[i]);
  }
  return out;
}

// RDKit's coordMap-tightened bounds matrix (setupInitialBoundsMatrix's success path with a coordMap):
// the raw topology bounds, then adjustBoundsMatFromCoordMap (pins each coordMap pair to its exact
// distance), then triangleSmoothBounds at tol 0.05. Flattened n*n row-major (upper above the
// diagonal, lower below), matching raw_bounds / native spec.bounds. Oracle for bb_embed::bounds::
// coord_map_bounds. `pin_idx[k]` is the atom index of pin k; `pin_xyz` is 3 doubles per pin. Empty
// on parse failure, size mismatch, or a triangle-smoothing failure (native's own fallback path,
// separately handled — see the test).
// RDKit's Stage-A force-field ENERGY (calcEnergy) at a supplied 4D geometry, with the
// firstMinimization FF (chiral 1.0, fourth 0.1, basinThresh). This is the exact value
// firstMinimization thresholds for its per-atom energy reject (MAX_MINIMIZED_E_PER_ATOM). NaN on
// failure. NB: RDKit's Chiral/FourthDim contribs use an energy inconsistent with their own gradient,
// so this can differ from native's gradient-consistent energy — the reason this oracle exists.
double stage_a_energy(rust::Str smiles, rust::Slice<const double> coords, double weight_chiral,
                      double weight_fourth, double basin) {
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) { return std::nan(""); }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) { RDKit::MolOps::findSSSR(*mol); }
  unsigned int n = mol->getNumAtoms();
  if (coords.size() != static_cast<std::size_t>(n) * 4) { return std::nan(""); }
  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  DistGeom::triangleSmoothBounds(mat);
  DistGeom::VECT_CHIRALSET chiralCenters, tetrahedralCenters;
  RDKit::DGeomHelpers::EmbeddingOps::findChiralSets(*mol, chiralCenters, tetrahedralCenters, nullptr);
  RDGeom::PointPtrVect positions;
  positions.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    auto *pt = new RDGeom::PointND(4);
    for (int c = 0; c < 4; ++c) { (*pt)[c] = coords[static_cast<std::size_t>(i) * 4 + c]; }
    positions.push_back(pt);
  }
  double e;
  {
    std::unique_ptr<ForceFields::ForceField> field(DistGeom::constructForceField(
        *mat, positions, chiralCenters, weight_chiral, weight_fourth, nullptr, basin, nullptr));
    field->initialize();
    e = field->calcEnergy();
  }
  for (auto *pt : positions) { delete pt; }
  return e;
}

rust::Vec<double> coord_map_bounds(rust::Str smiles, rust::Slice<const int> pin_idx,
                                   rust::Slice<const double> pin_xyz) {
  rust::Vec<double> out;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) { return out; }
  RDKit::MolOps::addHs(*mol);
  if (!mol->getRingInfo()->isInitialized()) { RDKit::MolOps::findSSSR(*mol); }
  std::size_t n = mol->getNumAtoms();
  if (pin_xyz.size() != pin_idx.size() * 3) { return out; }

  DistGeom::BoundsMatPtr mat = compute_raw_bounds(*mol);
  std::map<int, RDGeom::Point3D> coordMap;
  for (std::size_t k = 0; k < pin_idx.size(); ++k) {
    coordMap[pin_idx[k]] = RDGeom::Point3D(pin_xyz[k * 3 + 0], pin_xyz[k * 3 + 1], pin_xyz[k * 3 + 2]);
  }
  RDKit::DGeomHelpers::EmbeddingOps::adjustBoundsMatFromCoordMap(mat, (unsigned int)n, &coordMap);
  if (!DistGeom::triangleSmoothBounds(mat, 0.05)) {
    return out;  // smoothing failed — RDKit recomputes relaxed; native falls back to base (see test)
  }
  const double *data = mat->getData();
  out.reserve(n * n);
  for (std::size_t i = 0; i < n * n; ++i) { out.push_back(data[i]); }
  return out;
}

rust::Vec<double> smooth_bounds(rust::Slice<const double> bounds, std::size_t n, double tol) {
  DistGeom::BoundsMatrix mmat(static_cast<unsigned int>(n));
  double *d = mmat.getData();
  for (std::size_t i = 0; i < n * n; ++i) {
    d[i] = bounds[i];
  }
  DistGeom::triangleSmoothBounds(&mmat, tol);
  rust::Vec<double> out;
  out.reserve(n * n);
  for (std::size_t i = 0; i < n * n; ++i) {
    out.push_back(d[i]);
  }
  return out;
}

}  // namespace bb
