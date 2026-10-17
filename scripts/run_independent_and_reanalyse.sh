#!/usr/bin/env bash
# Independent section training followed by canonical minimum-change reevaluation.
# Requires verified external DLPFC .h5ad files and a working scientific environment.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OUT="${OUT:-${PWD}/LGA_FULL9_OUTPUT}"
export DATA_DIR="${DATA_DIR:?Set DATA_DIR to directory containing dlpfc_<section>.h5ad}"
export BASE="${BASE:-${PWD}}"
export BIB_FULL9_ROOT="$OUT"
export BIB_CANONICAL_OUT="${CANONICAL_OUT:-$OUT/canonical}"
mkdir -p "$OUT" "$BIB_CANONICAL_OUT"
bash "$REPO_ROOT/scripts/run_full_9pairs_3seeds.sh"
python3 "$REPO_ROOT/src/current/recompute_canonical.py"
python3 "$REPO_ROOT/src/current/verify_canonical_math.py"
printf '\nResults saved in %s and %s\n' "$OUT" "$BIB_CANONICAL_OUT"
