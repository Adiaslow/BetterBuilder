#!/bin/bash
# Build Boost into $BOOST_PREFIX using the toolchain's own GCC, so no C++ ABI boundary exists
# between Boost and the RDKit that links it. Idempotent.
set -euo pipefail
cd "$(dirname "$0")" && . ./versions.env

if [ -f "$BOOST_PREFIX/include/boost/version.hpp" ]; then
  echo "boost $BOOST_VERSION already installed at $BOOST_PREFIX"
  exit 0
fi

GXX="$GCC_PREFIX/bin/g++-${GCC_VERSION%.*}"
[ -x "$GXX" ] || { echo "build-gcc.sh must run first ($GXX missing)" >&2; exit 1; }

UNDERSCORED="${BOOST_VERSION//./_}"
SRC="$BB_BUILD_SCRATCH/boost"
mkdir -p "$SRC" && cd "$SRC"

[ -f "boost_$UNDERSCORED.tar.bz2" ] || \
  curl -fsSLO "https://archives.boost.io/release/$BOOST_VERSION/source/boost_$UNDERSCORED.tar.bz2"
[ -d "boost_$UNDERSCORED" ] || tar xf "boost_$UNDERSCORED.tar.bz2"
cd "boost_$UNDERSCORED"

# b2's config goes in scratch, not $HOME/user-config.jam — nothing lands outside the prefix.
# The rpath makes Boost's own libraries find the toolchain libstdc++ without LD_LIBRARY_PATH.
PYVER=$("$BB_PYTHON" -c "import sys;print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYINC="$BB_PYTHON_INCLUDE"
NPINC=$("$BB_PYTHON" -c "import numpy;print(numpy.get_include())")
cat > "$SRC/user-config.jam" <<EOF
using gcc : ${GCC_VERSION%.*} : $GXX : <linkflags>-Wl,-rpath,$GCC_PREFIX/lib64 ;
using python : $PYVER : $BB_PYTHON : $PYINC $NPINC ;
EOF

[ -x ./b2 ] || ./bootstrap.sh --prefix="$BOOST_PREFIX" --with-toolset=gcc

# serialization and iostreams are what RDKit links; the rest are what its CMake probes for.
./b2 -j"$BB_BUILD_JOBS" \
  --user-config="$SRC/user-config.jam" \
  --prefix="$BOOST_PREFIX" \
  toolset="gcc-${GCC_VERSION%.*}" \
  link=shared runtime-link=shared threading=multi variant=release \
  --with-serialization --with-iostreams --with-system --with-filesystem \
  --with-regex --with-thread --with-date_time --with-chrono --with-atomic \
  --with-program_options --with-python \
  install

grep -m1 "define BOOST_LIB_VERSION" "$BOOST_PREFIX/include/boost/version.hpp"
