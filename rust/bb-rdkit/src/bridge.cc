#include "bb-rdkit/src/bridge.h"

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
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/Atom.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

// Forward-declare the oracle's exported findDoubleBonds (in EmbeddingOps; no public header) so we
// call the real routine rather than reimplementing it. Symbol lives in libRDKitDistGeomHelpers.
namespace RDKit {
namespace DGeomHelpers {
namespace EmbeddingOps {
void findDoubleBonds(
    const RDKit::ROMol &mol,
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &doubleBondEnds,
    std::vector<std::pair<std::vector<unsigned int>, int>> &stereoDoubleBonds,
    const std::map<int, RDGeom::Point3D> *coordMap);
}  // namespace EmbeddingOps
}  // namespace DGeomHelpers
}  // namespace RDKit

namespace bb {

std::size_t num_heavy_atoms(rust::Str smiles) {
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  return mol ? static_cast<std::size_t>(mol->getNumHeavyAtoms()) : 0;
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

  // --- bounds matrix: exactly as the Python GetMoleculeBoundsMatrix (doTriangleSmoothing,
  //     useMacrocycle14config) ---
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
  RDKit::DGeomHelpers::initBoundsMat(mat);
  RDKit::DGeomHelpers::setTopolBounds(*mol, mat, /*set15bounds=*/true,
                                      /*scaleVDW=*/false,
                                      /*useMacrocycle14config=*/true);
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

  // --- chiral sets: reproduce Embedder.cpp findChiralSets (internal, not exported) using RDKit's
  //     own perception + bond-iteration order, so results are identical. vol (0,0) == tetrahedral
  //     (checks only). No coordMap here (base spec); atropisomers not yet handled. ---
  const RDKit::RingInfo *ri = mol->getRingInfo();
  for (const auto &atom : mol->atoms()) {
    if (atom->getAtomicNum() == 1) {
      continue;
    }
    auto ct = atom->getChiralTag();
    bool tagged = (ct == RDKit::Atom::CHI_TETRAHEDRAL_CW ||
                   ct == RDKit::Atom::CHI_TETRAHEDRAL_CCW);
    bool cn_deg4 = ((atom->getAtomicNum() == 6 || atom->getAtomicNum() == 7) &&
                    atom->getDegree() == 4);
    if (!(tagged || cn_deg4)) {
      continue;
    }
    std::vector<std::uint32_t> nbrs;
    RDKit::ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol->getAtomBonds(atom);
    while (beg != end) {
      nbrs.push_back((*mol)[*beg]->getOtherAtom(atom)->getIdx());
      ++beg;
    }
    if (nbrs.size() < 3) {
      continue;
    }
    double volLo = 5.0, volHi = 100.0;
    if (nbrs.size() < 4) {
      volLo = 2.0;
      nbrs.push_back(atom->getIdx());
    }
    unsigned int numSmallRings = 0;
    for (auto sz : ri->atomRingSizes(atom->getIdx())) {
      if (sz < 5) {
        ++numSmallRings;
      }
    }
    std::uint8_t fused = (numSmallRings > 1) ? 1 : 0;
    auto push = [&](double lo, double hi) {
      out->chi_center_.push_back(atom->getIdx());
      for (int k = 0; k < 4; ++k) {
        out->chi_atoms_.push_back(nbrs[k]);
      }
      out->chi_vol_.push_back(lo);
      out->chi_vol_.push_back(hi);
      out->chi_fused_.push_back(fused);
    };
    if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CCW) {
      push(volLo, volHi);
    } else if (ct == RDKit::Atom::CHI_TETRAHEDRAL_CW) {
      push(-volHi, -volLo);
    } else {  // untagged C/N deg-4 -> tetrahedral (checks only), with findChiralSets skip conditions
      if (ri->isInitialized() &&
          (ri->numAtomRings(atom->getIdx()) < 2 ||
           ri->isAtomInRingOfSize(atom->getIdx(), 3))) {
        // skip: <2 rings, or in a 3-ring
      } else {
        push(0.0, 0.0);
      }
    }
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
