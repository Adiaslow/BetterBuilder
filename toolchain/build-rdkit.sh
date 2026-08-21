#!/bin/bash
# Build the patched RDKit into $RDKIT_PREFIX: clone the pin, apply the patches from ../rdkit-patch,
# build against the toolchain's GCC, Boost and Eigen.
#
# Installs in-tree (RDKit's default), giving <prefix>/Code + <prefix>/lib — the layout
# BB_RDKIT_ROOT expects.
set -euo pipefail
cd "$(dirname "$0")" && . ./versions.env
PATCHES="$(cd .. && pwd)/rdkit-patch"

GXX="$GCC_PREFIX/bin/g++-${GCC_VERSION%.*}"
GCCBIN="$GCC_PREFIX/bin/gcc-${GCC_VERSION%.*}"
[ -x "$GXX" ] || { echo "build-gcc.sh must run first ($GXX missing)" >&2; exit 1; }
[ -f "$BOOST_PREFIX/include/boost/version.hpp" ] || { echo "build-boost.sh must run first" >&2; exit 1; }
[ -d "$EIGEN_PREFIX/include/eigen3" ] || { echo "build-eigen.sh must run first" >&2; exit 1; }

if [ ! -d "$RDKIT_PREFIX/.git" ]; then
  git clone https://github.com/rdkit/rdkit "$RDKIT_PREFIX"
fi
git -C "$RDKIT_PREFIX" checkout --quiet --force "$RDKIT_PIN"
git -C "$RDKIT_PREFIX" clean -qfd -e lib -e rdkit || true

git -C "$RDKIT_PREFIX" apply "$PATCHES/divya_amide_torsions.diff"
git -C "$RDKIT_PREFIX" apply "$PATCHES/cxx_compat.diff"
echo "patched:"; git -C "$RDKIT_PREFIX" diff --name-only | sed 's/^/  /'

NUMPY_INCLUDE=$("$BB_PYTHON" -c "import numpy;print(numpy.get_include())")

BUILD="$BB_BUILD_SCRATCH/rdkit"
rm -rf "$BUILD"; mkdir -p "$BUILD"; cd "$BUILD"

# Boost is found in CMake config mode, so it needs CMAKE_PREFIX_PATH — BOOST_ROOT alone drives
# the older module mode and is ignored here. RDKit requires Boost >= RDK_BOOST_VERSION (1.81).
cmake "$RDKIT_PREFIX" \
  -DCMAKE_C_COMPILER="$GCCBIN" \
  -DCMAKE_CXX_COMPILER="$GXX" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-mpopcnt" \
  -DCMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,$GCC_PREFIX/lib64" \
  -DCMAKE_SHARED_LINKER_FLAGS="-Wl,-rpath,$GCC_PREFIX/lib64 -Wl,-rpath,$BOOST_PREFIX/lib" \
  -DCMAKE_INSTALL_RPATH="$RDKIT_PREFIX/lib:$GCC_PREFIX/lib64:$BOOST_PREFIX/lib" \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -DCMAKE_PREFIX_PATH="$BOOST_PREFIX" \
  -DEIGEN3_INCLUDE_DIR="$EIGEN_PREFIX/include/eigen3" \
  -DRDK_BUILD_PYTHON_WRAPPERS=ON \
  -DPython3_EXECUTABLE="$BB_PYTHON" \
  -DPython3_INCLUDE_DIR="$BB_PYTHON_INCLUDE" \
  -DPython3_NumPy_INCLUDE_DIR="$NUMPY_INCLUDE" \
  -DPython3_FIND_STRATEGY=LOCATION \
  -DRDK_BUILD_CPP_TESTS=OFF \
  -DRDK_INSTALL_INTREE=ON \
  -DRDK_BUILD_COORDGEN_SUPPORT=ON \
  -DRDK_INSTALL_STATIC_LIBS=OFF

make -j"$BB_BUILD_JOBS"
make install

echo "installed $(ls "$RDKIT_PREFIX"/lib/*.so | wc -l) libraries in $RDKIT_PREFIX/lib"
