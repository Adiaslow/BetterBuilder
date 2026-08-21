#!/bin/sh
# Bench each pure-Rust component against C++ RDKit doing the same work, per component.
#
#   validation/bench/run.sh [corpus.smi]
#
# For every component we build to replace an RDKit piece, C++ RDKit is the speed floor: the Rust
# path has to be faster. Each C++ baseline calls RDKit's own routine (getExperimentalTorsions,
# setTopolBounds) — a stopwatch around the real library, not a reimplementation.
#
# Every harness reads the corpus up front, warms up, then times only the component; process startup
# and I/O are excluded on both sides so the ms/mol numbers compare like for like. Each prints an
# equal-work quantity (torsion count; constrained_pairs + checksum) that MUST agree across the Rust
# and C++ sides, or they are not doing the same work and the timing is meaningless.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
corpus=${1:-$repo/validation/seeds_5000.smi}

. "$repo/toolchain/env.sh"

R=$BB_RDKIT_ROOT
RDL="$R/External/RingFamilies/RingDecomposerLib/src/RingDecomposerLib"
workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

# Compile a C++ bench harness against the patched RDKit. One link line (superset of libraries) serves
# every harness.
build_cpp() {
  src=$1
  out=$2
  "$CXX" -O3 -std=c++20 "$src" -o "$out" \
    -I"$R/Code" -I"$RDL" -I"$BB_BOOST_PREFIX/include" -I"$BB_EIGEN_PREFIX/include/eigen3" \
    -L"$R/lib" -L"$BB_BOOST_PREFIX/lib" \
    -lRDKitSmilesParse -lRDKitGraphMol -lRDKitRDGeneral -lRDKitRDGeometryLib \
    -lRDKitDistGeomHelpers -lRDKitDistGeometry \
    -lRDKitForceFieldHelpers -lRDKitForceField -lRDKitSubstructMatch -lRDKitDataStructs \
    -Wl,-rpath,"$R/lib" -Wl,-rpath,"$BB_BOOST_PREFIX/lib" -Wl,-rpath,"$BB_CXX_PREFIX/lib64"
}

# Build a Rust bench binary if it is declared; return non-zero (and leave a note) if not yet built.
rust_bin() {
  pkg=$1
  name=$2
  if cargo build --release -q -p "$pkg" --bin "$name" \
      --manifest-path "$repo/rust/Cargo.toml" 2>/dev/null; then
    echo "$repo/rust/target/release/$name"
    return 0
  fi
  return 1
}

echo "building C++ baselines against patched RDKit at $R ..."
build_cpp "$here/time_cpp_torsions.cc" "$workdir/time_cpp_torsions"
build_cpp "$here/time_cpp_bounds.cc" "$workdir/time_cpp_bounds"
build_cpp "$here/time_cpp_uff.cc" "$workdir/time_cpp_uff"

echo
echo "corpus: $corpus ($(grep -c . "$corpus") lines)"

echo
echo "=== UFF atom typing ==="
if rbin=$(rust_bin bb-bounds bench-uff); then
  echo "--- Rust (bb-bounds) ---"
  "$rbin" "$corpus"
fi
echo "--- C++ RDKit (floor) ---"
"$workdir/time_cpp_uff" "$corpus"

echo
echo "=== torsion assignment ==="
if rbin=$(rust_bin bb-perceive bench-torsions); then
  echo "--- Rust (bb-perceive) ---"
  "$rbin" "$corpus"
fi
echo "--- C++ RDKit (floor) ---"
"$workdir/time_cpp_torsions" "$corpus"

echo
echo "=== setTopolBounds ==="
if rbin=$(rust_bin bb-bounds bench-bounds); then
  echo "--- Rust (bb-bounds) ---"
  "$rbin" "$corpus"
else
  echo "--- Rust (bb-bounds): pending — constructor not yet built (milestones M1-M5) ---"
fi
echo "--- C++ RDKit (floor) ---"
"$workdir/time_cpp_bounds" "$corpus"
