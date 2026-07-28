#!/usr/bin/env bash
# Assert that TempLat's lattice kernels compile to VECTOR code on the CPU.
#
# Why this exists: a scalar kernel is invisible in every other test -- it computes exactly the right
# answer, just ~2.5x too slowly -- which is how TempLat shipped a 100%-scalar SU(2) kick for as long
# as it did. See benchmarks/PERFORMANCE.md. The whole fix is one independence assertion on the
# innermost loop (TEMPLAT_ASSUME_INDEPENDENT, in parallel/devices/kokkos/kokkos.h); if that loop ever
# stops being visible to the vectorizer -- a refactor of the dispatch, a data-dependent branch in the
# kernel body, a compiler that stops honouring the pragma -- the packed instructions silently
# disappear and only this check notices.
#
# Usage: check_vectorization.sh [--cpu-build] <binary> [min_packed]
# Skips (exit 0) on anything that isn't a host x86-64 build with objdump available.
#
# --cpu-build says the caller knows this binary was compiled for the CPU, which makes the ABSENCE of
# the CPU kernels a failure rather than a skip. Without it, a revert of the dispatch would look like a
# GPU build and quietly pass -- which is the exact failure this check is here to catch. CTest passes
# it; a hand-run against some arbitrary (possibly GPU) binary does not.

set -uo pipefail

CPU_BUILD=0
if [ "${1:-}" = "--cpu-build" ]; then
  CPU_BUILD=1
  shift
fi

BIN=${1:?usage: check_vectorization.sh [--cpu-build] <binary> [min_packed]}
MIN_PACKED=${2:-100}

if ! command -v objdump >/dev/null 2>&1; then
  echo "SKIP: objdump not available"
  exit 0
fi
if [ "$(uname -m)" != "x86_64" ]; then
  echo "SKIP: packed-instruction patterns are x86-64 specific (found $(uname -m))"
  exit 0
fi
if [ ! -f "$BIN" ]; then
  echo "FAIL: no such binary: $BIN"
  exit 1
fi

# The CPU lattice dispatch runs every field assignment through KokkosNDLambdaWrapperInnerLoop, so the
# hot kernels are exactly the symbols that mention it. Count packed (…pd) vs scalar (…sd) FP ops in
# those, and nowhere else -- a binary-wide count would be satisfied by some other loop vectorizing.
read -r PACKED SCALAR < <(
  objdump -d --no-show-raw-insn "$BIN" | awk '
    /^[0-9a-f]+ <.*>:/ { hot = ($0 ~ /KokkosNDLambdaWrapperInnerLoop/) }
    hot && $1 ~ /^[0-9a-f]+:$/ {
      if ($2 ~ /^v(fmadd|fmsub|fnmadd|fnmsub|mul|add|sub|div)[0-9]*pd$/) p++
      else if ($2 ~ /^v?(fmadd|fmsub|fnmadd|fnmsub|mul|add|sub|div)[0-9]*sd$/) s++
    }
    END { print p+0, s+0 }'
)

if [ "$PACKED" -eq 0 ] && [ "$SCALAR" -eq 0 ]; then
  if [ "$CPU_BUILD" -eq 0 ]; then
    echo "SKIP: no KokkosNDLambdaWrapperInnerLoop kernels in $(basename "$BIN") (GPU build, or nothing inlined)"
    exit 0
  fi
  cat <<EOF
FAIL: $(basename "$BIN") is a CPU build, but contains no KokkosNDLambdaWrapperInnerLoop kernels at all.

The CPU lattice dispatch is gone: device::iteration::foreach is no longer routing field assignments
through the inner-loop wrapper, so nothing can vectorize. See
include/TempLat/parallel/devices/kokkos/kokkos_iteration.h and benchmarks/PERFORMANCE.md.
EOF
  exit 1
fi

echo "$(basename "$BIN"): packed=$PACKED scalar=$SCALAR (in the CPU lattice kernels)"

if [ "$PACKED" -lt "$MIN_PACKED" ]; then
  cat <<EOF
FAIL: the lattice kernels are (almost) all scalar -- $PACKED packed FP ops, expected >= $MIN_PACKED.

The vectorizer has stopped seeing the innermost lattice loop. It will not tell you this itself: every
correctness test still passes, the code is just ~2.5x slower. Start at TEMPLAT_ASSUME_INDEPENDENT and
KokkosNDLambdaWrapperInnerLoop (include/TempLat/parallel/devices/kokkos/), and read
benchmarks/PERFORMANCE.md. Reproduce with:
  g++ -O3 -march=native -fopt-info-vec-missed ...
EOF
  exit 1
fi

echo "PASS"
