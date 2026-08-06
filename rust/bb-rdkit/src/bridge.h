#pragma once
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "rust/cxx.h"

namespace bb {

// Holds the extracted distance-geometry problem for one molecule. Opaque to Rust;
// accessed via the const methods below (cxx maps std::vector<T> <-> CxxVector<T>).
class MolSpecC {
 public:
  bool ok_ = false;
  std::size_t n_atoms_ = 0;
  std::vector<double> bounds_;            // n_atoms*n_atoms row-major; upper-tri=UB, lower-tri=LB
  std::vector<std::uint32_t> tors_atoms_; // 4 per experimental torsion
  std::vector<double> tors_v_;            // 6 per torsion (force constants)
  std::vector<std::int32_t> tors_signs_;  // 6 per torsion
  // chiral sets (reproduced from RDKit findChiralSets). vol (0,0) => tetrahedral, checks only.
  std::vector<std::uint32_t> chi_center_; // 1 per set (center atom / centroid idx)
  std::vector<std::uint32_t> chi_atoms_;  // 4 per set (the volume points)
  std::vector<double> chi_vol_;           // 2 per set (lower, upper)
  std::vector<std::uint8_t> chi_fused_;   // 1 per set (in fused small rings)
  std::vector<std::uint32_t> imp_atoms_;  // 4 per UFF improper contrib (i, j=center, k, l)
  std::vector<double> imp_coef_;          // 4 per contrib (C0, C1, C2, fc [incl. ×10 scaling])
  std::vector<std::uint32_t> bonds_;      // 2 per bond (begin, end)
  std::vector<std::uint32_t> angles_;     // 4 per angle (i, j=center, k, tripleFlag)
  std::vector<std::uint32_t> pin_atoms_;  // core-pin recipe: atoms frozen to the seed conformer
  std::uint32_t core_seeds_ = 0;          // rdkit_confs_1 (20 if exo<=2 else 10)
  std::uint32_t sidechain_confs_ = 0;     // per-seed confs (10 if exo<=2 else 20)
  std::vector<std::uint32_t> dbe_;        // 3 per double-bond end (nbr, dbAtom, otherDbAtom)
  std::vector<std::uint32_t> sdb_atoms_;  // 4 per stereo double bond (s0, begin, end, s1)
  std::vector<std::int32_t> sdb_sign_;    // 1 per stereo double bond (+1 trans/E, -1 cis/Z)

  bool ok() const { return ok_; }
  std::size_t n_atoms() const { return n_atoms_; }
  const std::vector<double> &bounds() const { return bounds_; }
  const std::vector<std::uint32_t> &tors_atoms() const { return tors_atoms_; }
  const std::vector<double> &tors_v() const { return tors_v_; }
  const std::vector<std::int32_t> &tors_signs() const { return tors_signs_; }
  const std::vector<std::uint32_t> &chi_center() const { return chi_center_; }
  const std::vector<std::uint32_t> &chi_atoms() const { return chi_atoms_; }
  const std::vector<double> &chi_vol() const { return chi_vol_; }
  const std::vector<std::uint8_t> &chi_fused() const { return chi_fused_; }
  const std::vector<std::uint32_t> &imp_atoms() const { return imp_atoms_; }
  const std::vector<double> &imp_coef() const { return imp_coef_; }
  const std::vector<std::uint32_t> &bonds() const { return bonds_; }
  const std::vector<std::uint32_t> &angles() const { return angles_; }
  const std::vector<std::uint32_t> &pin_atoms() const { return pin_atoms_; }
  std::uint32_t core_seeds() const { return core_seeds_; }
  std::uint32_t sidechain_confs() const { return sidechain_confs_; }
  const std::vector<std::uint32_t> &dbe() const { return dbe_; }
  const std::vector<std::uint32_t> &sdb_atoms() const { return sdb_atoms_; }
  const std::vector<std::int32_t> &sdb_sign() const { return sdb_sign_; }
};

std::size_t num_heavy_atoms(rust::Str smiles);
std::unique_ptr<MolSpecC> build_mol_spec(rust::Str smiles);

// Test-only: run RDKit's real triangleSmoothBounds on a flat n×n bounds matrix (row<col = upper,
// row>col = lower — same layout as bb-core). Used to validate bb-embed's Rust port bit-for-bit.
rust::Vec<double> smooth_bounds(rust::Slice<const double> bounds, std::size_t n, double tol);

}  // namespace bb
