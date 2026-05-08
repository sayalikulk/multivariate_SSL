#!/usr/bin/env bash
# Build & install the master and quic-improvements versions of mSSL into
# separate library trees, then run the R benchmark that loads each in turn.
#
# Usage:   bash tests/run_benchmark.sh
# Outputs: tests/results/{base,improved}/ + tests/results/comparison.txt

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORKTREE="$ROOT/.bench-master"
LIB_BASE="$ROOT/tests/lib_base"
LIB_NEW="$ROOT/tests/lib_improved"
RESULTS="$ROOT/tests/results"

mkdir -p "$LIB_BASE" "$LIB_NEW" "$RESULTS"

echo "==> [1/4] Preparing master worktree at $WORKTREE"
if [ ! -d "$WORKTREE" ]; then
  git -C "$ROOT" worktree add "$WORKTREE" master
else
  git -C "$WORKTREE" fetch origin master --quiet || true
  git -C "$WORKTREE" checkout master --quiet
fi

echo "==> [2/4] Installing BASE (master) into $LIB_BASE"
R CMD INSTALL --no-multiarch --library="$LIB_BASE" "$WORKTREE/mSSL"

echo "==> [3/4] Installing IMPROVED (current branch) into $LIB_NEW"
R CMD INSTALL --no-multiarch --library="$LIB_NEW" "$ROOT/mSSL"

echo "==> [4/4] Running benchmark"
Rscript "$ROOT/tests/benchmark.R" "$LIB_BASE" "$LIB_NEW" "$RESULTS"

echo
echo "Done. See $RESULTS/comparison.txt"
