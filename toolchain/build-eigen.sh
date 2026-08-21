#!/bin/bash
# Install Eigen headers into $EIGEN_PREFIX. Header-only, so this unpacks rather than compiles.
set -euo pipefail
cd "$(dirname "$0")" && . ./versions.env

if [ -d "$EIGEN_PREFIX/include/eigen3/Eigen" ]; then
  echo "eigen $EIGEN_VERSION already installed at $EIGEN_PREFIX"
  exit 0
fi

SRC="$BB_BUILD_SCRATCH/eigen"
mkdir -p "$SRC" "$EIGEN_PREFIX/include/eigen3" && cd "$SRC"

[ -f "eigen-$EIGEN_VERSION.tar.gz" ] || \
  curl -fsSLO "https://gitlab.com/libeigen/eigen/-/archive/$EIGEN_VERSION/eigen-$EIGEN_VERSION.tar.gz"
[ -d "eigen-$EIGEN_VERSION" ] || tar xf "eigen-$EIGEN_VERSION.tar.gz"

cp -r "eigen-$EIGEN_VERSION/Eigen" "eigen-$EIGEN_VERSION/unsupported" "$EIGEN_PREFIX/include/eigen3/"
echo "eigen $EIGEN_VERSION installed at $EIGEN_PREFIX"
