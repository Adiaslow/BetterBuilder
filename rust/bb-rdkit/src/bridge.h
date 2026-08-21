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
  std::vector<std::uint8_t> atomic_nums_; // 1 per atom, in addHs order
  std::int32_t formal_charge_ = 0;        // net formal charge of the molecule
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
  const std::vector<std::uint8_t> &atomic_nums() const { return atomic_nums_; }
  std::int32_t formal_charge() const { return formal_charge_; }
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

// Oracle for the pure-Rust setTopolBounds port: RDKit's bounds matrix straight out of
// initBoundsMat + setTopolBounds, BEFORE triangleSmoothBounds. The production bridge only exposes
// the smoothed matrix, so this exposes the raw one to gate the constructor in isolation. Same
// molecule construction and same setTopolBounds flags as build_mol_spec (both route through
// compute_raw_bounds), so smoothing this reproduces build_mol_spec's bounds exactly. Returns a flat
// n*n row-major matrix (upper-tri = UB, lower-tri = LB), or empty on a SMILES RDKit cannot parse.
rust::Vec<double> raw_bounds(rust::Str smiles);

// RDKit's UFF parameter table (ForceFields::UFF::defaultParamData) verbatim, for vendoring. This is
// the exact string setTopolBounds' UFF atom typing parses, so the Rust port reads a byte-identical
// copy rather than a transcription.
rust::String uff_param_table();

// Per-atom UFF atom-type labels (RDKit::UFF::Tools::getAtomLabel) for the molecule, in post-AddHs
// order — the oracle for the Rust getAtomLabel port. Empty on a SMILES RDKit cannot parse.
rust::Vec<rust::String> uff_atom_labels(rust::Str smiles);

// Per-bond UFF rest length exactly as set12Bounds computes accumData.bondLengths (calcBondRestLength
// from the two atoms' UFF params, or (vdw1+vdw2)/2 when either lacks params), in post-AddHs bond
// order — the oracle for the Rust calcBondRestLength port. Empty on a SMILES RDKit cannot parse.
rust::Vec<double> uff_bond_rest_lengths(rust::Str smiles);

// RDKit's topological distance matrix (MolOps::getDistanceMat, useBO=false, useAtomWts=false), flat
// n*n row-major over post-AddHs atoms — the oracle for the Rust distance_matrix port. Empty on a
// SMILES RDKit cannot parse.
rust::Vec<double> topo_distance_matrix(rust::Str smiles);

// RDKit's per-bond double-bond stereo as SmilesToMol assigned it: three ints per bond in post-AddHs
// bond order — getStereo() enum, then the two getStereoAtoms (-1 when unset). The oracle for the
// Rust double-bond-stereo port. Empty on a SMILES RDKit cannot parse.
rust::Vec<int> set_arom_only(rust::Str smiles);
rust::Vec<int> atom_elec(rust::Str smiles);
rust::Vec<int> aromatic_perception(rust::Str smiles);
rust::Vec<int> bond_dirs(rust::Str smiles);
rust::Vec<double> stage_a_ff_grad(rust::Str smiles, rust::Slice<const double> coords, double weight_chiral, double weight_fourth, double basin);
rust::Vec<double> stage_c_ff_grad(rust::Str smiles, rust::Slice<const double> coords);
rust::Vec<double> stage_c_improper_grad(rust::Str smiles, rust::Slice<const double> coords);
rust::Vec<double> stage_c_torsion_grad(rust::Str smiles, rust::Slice<const double> coords);
// RDKit's post-embed accept/reject checks evaluated on the given 3D conformer, in native
// passes_checks order: 6 ints {tetrahedral, chiral, planarity, double-bond-geometry, final-chiral,
// double-bond-stereo}, each 1=pass 0=fail. Empty on parse/size failure.
rust::Vec<int> embed_checks(rust::Str smiles, rust::Slice<const double> coords);
// RDKit's coordMap-tightened bounds matrix (setupInitialBoundsMatrix success path): raw topology
// bounds + adjustBoundsMatFromCoordMap + triangleSmoothBounds(tol=0.05), flattened n*n row-major.
// Empty on parse failure, size mismatch, or smoothing failure. Oracle for coord_map_bounds.
rust::Vec<double> coord_map_bounds(rust::Str smiles, rust::Slice<const int> pin_idx,
                                   rust::Slice<const double> pin_xyz);
// RDKit's powerEigenSolver on a caller-supplied metric matrix T (lower triangle incl. diagonal,
// size n(n+1)/2). Returns [n_eig eigenvalues][n_eig*n eigenvector components]; empty on size
// mismatch. Isolation oracle for native's exact metric-matrix eigenstep. See bridge.cc.
// RDKit's Stage-A calcEnergy at a 4D geometry (firstMinimization FF). The value the per-atom energy
// reject thresholds. NaN on failure. See bridge.cc.
double stage_a_energy(rust::Str smiles, rust::Slice<const double> coords, double weight_chiral,
                      double weight_fourth, double basin);
rust::Vec<int> bond_types(rust::Str smiles);
rust::Vec<int> sssr_rings(rust::Str smiles);
rust::Vec<int> bond_stereo(rust::Str smiles);

// RDKit's per-atom getChiralTag() in post-AddHs order (CHI_UNSPECIFIED=0, CW=1, CCW=2, ...) — the
// oracle for the Rust tetrahedral-chirality port. Empty on a SMILES RDKit cannot parse.
rust::Vec<int> chiral_tags(rust::Str smiles);

// RDKit's raw (pre-smoothing) bounds matrix with set13/set14/set15 selectable, flat n*n row-major —
// the staged oracle for the bounds-matrix port (set13=set14=set15=false is set12+VDW only). Other
// flags match compute_raw_bounds. Empty on a SMILES RDKit cannot parse.
rust::Vec<double> raw_bounds_stage(rust::Str smiles, bool set13, bool set14, bool set15);

}  // namespace bb
