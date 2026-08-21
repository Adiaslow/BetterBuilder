#!/bin/bash
# Capture oracle stage-level state for differential testing.
#
#   BB_ORACLE_CONTAINER=$HOME/toolchain/images/build_macrocycle_final.sif \
#   BB_ORACLE_RUN=$HOME/toolchain/oracle-run/full/out \
#   validation/parity/dump.sh [count]
#
# BB_ORACLE_RUN is the working directory of a completed pipeline run; its solv/<mol>/output.{mol2,solv}
# are the oracle's own intermediates. Output lands in validation/parity/fixtures/<molecule>.json.
#
# Nothing here modifies the oracle: dump_oracle.py imports its modules and reads what they compute.
set -euo pipefail
cd "$(dirname "$0")/../.."
COUNT="${1:-5}"
OUT=validation/parity/fixtures

: "${BB_ORACLE_CONTAINER:?set BB_ORACLE_CONTAINER to the production image (build_macrocycle_final.sif)}"
: "${BB_ORACLE_RUN:?set BB_ORACLE_RUN to a completed pipeline run directory}"

mkdir -p "$OUT"
n=0
for d in "$BB_ORACLE_RUN"/solv/*/; do
  [ "$n" -ge "$COUNT" ] && break
  [ -f "$d/output.mol2" ] && [ -f "$d/output.solv" ] || continue
  id=$(basename "$d"); id="${id%%_*}"
  apptainer exec \
      -B "$BB_ORACLE_RUN":/run \
      -B "$(pwd)/validation/parity":/parity \
      "$BB_ORACLE_CONTAINER" \
      python3 /parity/dump_oracle.py \
          "/run/solv/$(basename "$d")/output.mol2" \
          "/run/solv/$(basename "$d")/output.solv" \
          "/parity/fixtures/$id.json" \
    | sed 's/^/  /'
  n=$((n+1))
done
echo "wrote $n oracle state dumps to $OUT"
