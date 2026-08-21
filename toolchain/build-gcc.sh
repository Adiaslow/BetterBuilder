#!/bin/bash
# Build GCC into $GCC_PREFIX. Idempotent: exits early if the compiler is already installed.
set -euo pipefail
cd "$(dirname "$0")" && . ./versions.env

if [ -x "$GCC_PREFIX/bin/g++-${GCC_VERSION%.*}" ]; then
  echo "gcc $GCC_VERSION already installed at $GCC_PREFIX"
  exit 0
fi

SUFFIX="-${GCC_VERSION%.*}"
SRC="$BB_BUILD_SCRATCH/gcc"
mkdir -p "$SRC" "$GCC_PREFIX"
cd "$SRC"

[ -f "gcc-$GCC_VERSION.tar.xz" ] || curl -fLO "https://ftp.gnu.org/gnu/gcc/gcc-$GCC_VERSION/gcc-$GCC_VERSION.tar.xz"
[ -d "gcc-$GCC_VERSION" ] || tar xf "gcc-$GCC_VERSION.tar.xz"
[ -d "gcc-$GCC_VERSION/gmp" ] || ( cd "gcc-$GCC_VERSION" && ./contrib/download_prerequisites )

mkdir -p build && cd build
[ -f Makefile ] || "../gcc-$GCC_VERSION/configure" \
  --prefix="$GCC_PREFIX" \
  --enable-languages=c,c++ \
  --disable-multilib \
  --with-system-zlib \
  --program-suffix="$SUFFIX"

# Full 3-stage bootstrap: stages 2 and 3 must produce identical compilers, which self-checks
# the build rather than trusting it.
make -j"$BB_BUILD_JOBS"
make install

"$GCC_PREFIX/bin/g++$SUFFIX" --version | head -1
