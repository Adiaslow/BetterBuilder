#!/bin/sh
# Reproducible entry point for the BetterBuilder gate suite.
#
#   rust/gates.sh                 # everything: fast + live + verify (corpus = seeds_5000)
#   rust/gates.sh fast            # ONLY toolchain-free gates (native vs committed goldens) — no RDKit
#   rust/gates.sh verify          # regenerate goldens from live RDKit, byte-diff committed (anti-stale)
#   rust/gates.sh certify         # prove our from-source build == her container build (byte-identical)
#   rust/gates.sh live seeds_100.smi   # live RDKit-parity gates on a chosen corpus
#   rust/gates.sh hunt            # divergence-hunt: generate adversarial macrocycles, diff native vs bridge
#   rust/gates.sh dbgate          # db2 writer corpus gate: native per-block db2 == her mol2db2 (byte-identical); full pipeline with AMSOL
#   rust/gates.sh bench           # speed gate: each pure-Rust component beats the C++ RDKit floor (equal work)
#
# Sources the committed toolchain env (toolchain/env.sh) — never a scratch copy. A red gate aborts
# non-zero: red is a signal to fix the code, never licence to weaken or delete a gate.
set -eu

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO/rust"

MODE="${1:-all}"
CORPUS_ARG="${2:-$REPO/validation/seeds_5000.smi}"
case "$CORPUS_ARG" in /*) CORPUS="$CORPUS_ARG" ;; *) CORPUS="$REPO/$CORPUS_ARG" ;; esac

# Oracle container — ONE definition (env-var BB_ORACLE_CONTAINER), shared by certify + hunt. The
# per-mode `[ -f "$IMG" ]` existence check stays in those functions so fast/verify don't require it.
APP="${APPTAINER:-/usr/bin/apptainer}"
IMG="${BB_ORACLE_CONTAINER:-/nfs/home/amurray2/toolchain/images/build_macrocycle_final.sif}"

# --- Fast tier: native correctness + parity vs committed goldens. No RDKit toolchain. ---
run_fast() {
  echo "== fast tier (no RDKit toolchain): native correctness + golden parity =="
  cargo test -p bb-embed --test embedding      -- --nocapture
  cargo test -p bb-embed --test parity_goldens -- --nocapture
}

# --- Verify tier: regenerate goldens from LIVE RDKit and byte-diff the committed ones. This is the
#     anti-staleness / version-bump guard; the committed goldens must equal a fresh live run. ---
run_verify() {
  echo "== verify tier: committed goldens == live RDKit ($(cat "$REPO/validation/fixtures/rdkit/VERSION")) =="
  . "$REPO/toolchain/env.sh"
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' EXIT
  cargo run -q -p bb-rdkit --bin gen-fixtures -- "$REPO/validation/seeds_100.smi" "$tmp"
  for f in VERSION stage_a_grad.jsonl; do
    if ! diff -q "$tmp/$f" "$REPO/validation/fixtures/rdkit/$f" >/dev/null; then
      echo "STALE: validation/fixtures/rdkit/$f differs from live RDKit — regenerate (see fixtures/README.md)" >&2
      exit 1
    fi
  done
  echo "goldens match live RDKit"
}

# --- Certify tier (CONFORMANCE.md Phase 1): prove our from-source RDKit build == Divya's actual
#     production build, by compiling oracle dumpers against HER container libs and diffing vs goldens.
#     Byte-identical or it fails. Needs her container (BB_ORACLE_CONTAINER) + the toolchain env. ---
run_certify() {
  echo "== certify tier: our from-source build == her container build (byte-identical) =="
  . "$REPO/toolchain/env.sh"
  [ -f "$IMG" ] || { echo "oracle container not found: $IMG (set BB_ORACLE_CONTAINER)" >&2; exit 2; }
  CDIR="$REPO/validation/certify"
  # The repo lives at $REPO, which apptainer does not auto-mount, so bind it into the container at its
  # own path (rw): the committed sources compile in place and the binary lands in a persistent,
  # gitignored build dir under the repo — no ephemeral /tmp. The toolchain headers are under $HOME,
  # which apptainer auto-mounts. `BB_ORACLE_BIND` carries the same bind to certify.py's run step.
  BIND="$REPO:$REPO"
  BUILD="$CDIR/.build"; mkdir -p "$BUILD"
  inc="-I$BB_RDKIT_ROOT/Code \
       -I$BB_RDKIT_ROOT/External/RingFamilies/RingDecomposerLib/src/RingDecomposerLib \
       -I$BB_BOOST_PREFIX/include -I$BB_EIGEN_PREFIX/include/eigen3"
  libs="-L/soft/rdkit_libs -Wl,-rpath,/soft/rdkit_libs \
        -lRDKitSmilesParse -lRDKitGraphMol -lRDKitRDGeneral -lRDKitRDGeometryLib \
        -lRDKitDistGeomHelpers -lRDKitDistGeometry -lRDKitForceFieldHelpers \
        -lRDKitForceField -lRDKitSubstructMatch"
  # compile each distinct dumper source once against her libs
  for stem in stage_a_dump stage_c_dump checks_dump coord_map_dump; do
    $APP exec --bind "$BIND" "$IMG" g++ -std=c++20 -O2 "$CDIR/$stem.cc" $inc $libs -o "$BUILD/$stem" \
      || { echo "compile FAILED: $stem against her libs" >&2; exit 1; }
  done
  # certify each oracle: <binary> <golden> [dumper args…]. Stage-A/B share stage_a_dump (constructForceField),
  # differing only in weights: A = firstMinimization (1.0 0.1), B = minimizeFourthDimension (0.2 1.0).
  # coord_map_dump serves both bounds (k=0 pins) and coord_map (embedded-conformer pins).
  certify() {
    bin="$1"; gold="$2"; shift 2
    BB_ORACLE_CONTAINER="$IMG" BB_ORACLE_BIND="$BIND" python3 "$CDIR/certify.py" 100 "$BUILD/$bin" \
      "$REPO/validation/fixtures/rdkit/$gold.jsonl" "$@" || exit 1
  }
  certify stage_a_dump    stage_a_grad 1.0 0.1
  certify stage_a_dump    stage_b_grad 0.2 1.0
  certify stage_c_dump    stage_c_grad
  certify checks_dump     checks_grad
  certify coord_map_dump  bounds_grad
  certify coord_map_dump  coord_map_grad
}

# --- Live tier: RDKit-parity gates not yet converted to goldens + the corpus-scale gate. ---
run_live() {
  echo "== live RDKit-parity gates ($CORPUS) =="
  . "$REPO/toolchain/env.sh"
  cargo test -p bb-rdkit --test ff_parity
  cargo test -p bb-rdkit --test embed_checks_parity
  cargo test -p bb-rdkit --test coord_map_bounds_parity
  BB_PARITY_CORPUS="$CORPUS" cargo test -p bb-rdkit --test corpus_parity -- --ignored --nocapture
}

# --- Hunt tier (CONFORMANCE.md Phase 2): divergence-hunt the deterministic layers. Generate adversarial
#     macrocycles off her real workload and diff native vs the RDKit bridge (spec layer) + run the
#     vetted FF/checks/bounds gate on them (embed layer). "no diffs" means "searched hard." ---
run_hunt() {
  echo "== hunt tier: divergence-hunting native vs RDKit bridge on generated adversarial macrocycles =="
  [ -f "$IMG" ] || { echo "oracle container not found: $IMG (set BB_ORACLE_CONTAINER)" >&2; exit 2; }
  H="$REPO/validation/hunt"; GEN="$H/generated.smi"
  NPER="${BB_HUNT_NPER:-2}"; HSEED="${BB_HUNT_SEED:-20210185}"
  $APP exec --bind "$REPO:$REPO" "$IMG" python3 "$H/mutate_corpus.py" "$CORPUS" "$GEN" "$NPER" "$HSEED" || exit 1
  APPTAINERENV_HUNT_NPROC="${HUNT_NPROC:-24}" \
    $APP exec --bind "$REPO:$REPO" "$IMG" python3 "$H/differential.py" "$GEN" || exit 1
  # embed layer: point the vetted corpus gate at a bounded subset of the generated corpus (host toolchain).
  head -"${BB_HUNT_EMBED_N:-500}" "$GEN" > "$H/generated_embed.smi"
  ( . "$REPO/toolchain/env.sh"
    BB_PARITY_CORPUS="$H/generated_embed.smi" cargo test --release -p bb-rdkit --test corpus_parity -- --ignored --nocapture ) || exit 1
}


# --- Output-writer corpus gate: native's per-block db2 byte-identical to her container mol2db2 (stage-
#     isolated, same mol2+solv → both writers). Needs AMSOL (host) + her container. ---
run_dbgate() {
  echo "== db2 writer corpus gate: native per-block db2 == her mol2db2 (byte-identical modulo SMILES metadata) =="
  python3 "$REPO/validation/db2/corpus_writer_gate.py" "$REPO/validation/seeds_100.smi" "${BB_DBGATE_N:-30}"
  # the full pipeline with real AMSOL charges (ignored by default — it needs AMSOL, as this tier does)
  cargo test -p bb-output --test full_pipeline -- --ignored
}

# --- Bench tier: the SPEED gate. `run.sh` times each pure-Rust component against C++ RDKit doing the
#     same work; `bench_gate.py` FAILS if the work quantities disagree (timing meaningless) or if Rust
#     is not under the C++ floor. C++ RDKit is the floor [[speed-floor-is-cpp-rdkit]]. Needs the
#     toolchain env (run.sh sources it) to compile the C++ baselines. ---
run_bench() {
  echo "== bench tier: each pure-Rust component beats the C++ RDKit floor doing equal work =="
  bcorpus="${BB_BENCH_CORPUS:-$REPO/validation/seeds_100.smi}"
  case "$bcorpus" in /*) ;; *) bcorpus="$REPO/$bcorpus" ;; esac  # resolve relative to the repo, like $CORPUS
  out="$(mktemp)"
  # Merge stderr: the Rust harnesses print their timing line to stderr, the C++ ones to stdout — the
  # gate needs both halves, so 2>&1 before the pipe (else it sees "MISSING rust half").
  sh "$REPO/validation/bench/run.sh" "$bcorpus" 2>&1 | tee "$out"
  python3 "$REPO/validation/bench/bench_gate.py" "$out" || { rm -f "$out"; exit 1; }
  rm -f "$out"
}

case "$MODE" in
  fast)    run_fast ;;
  verify)  run_verify ;;
  certify) run_certify ;;
  live)    run_live ;;
  hunt)    run_hunt ;;
  dbgate)  run_dbgate ;;
  bench)   run_bench ;;
  all)     run_fast; run_verify; run_certify; run_live; run_hunt; run_dbgate; run_bench ;;
  *) echo "usage: gates.sh [fast|verify|certify|live|hunt|dbgate|bench|all] [corpus.smi]" >&2; exit 2 ;;
esac
echo "GATES PASSED ($MODE)"
