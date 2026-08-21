# Point a cargo build at this toolchain:  . toolchain/env.sh
#
# Exports only what bb-rdkit's build script reads. Prefixes are per-package, so no combined
# directory has to be assembled.

. "$(dirname "${BASH_SOURCE[0]:-$0}")/versions.env"

export BB_RDKIT_ROOT="$RDKIT_PREFIX"
export BB_BOOST_PREFIX="$BOOST_PREFIX"
export BB_EIGEN_PREFIX="$EIGEN_PREFIX"
# The C++ runtime the RDKit libraries were built against. bb-rdkit's build script links and
# rpaths this, so the system libstdc++ (which may be older) is not used by accident.
export BB_CXX_PREFIX="$GCC_PREFIX"
export CXX="$GCC_PREFIX/bin/g++-${GCC_VERSION%.*}"
export CC="$GCC_PREFIX/bin/gcc-${GCC_VERSION%.*}"
