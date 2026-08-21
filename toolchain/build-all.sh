#!/bin/bash
# Build the whole toolchain in dependency order. Each step is idempotent, so re-running after a
# failure resumes rather than restarting.
set -euo pipefail
cd "$(dirname "$0")"

for step in gcc eigen boost rdkit; do
  echo "===== $step ====="
  ./build-$step.sh
done

echo "===== done ====="
echo "point a build at it with:  . toolchain/env.sh"
