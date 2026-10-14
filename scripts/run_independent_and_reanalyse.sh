#!/usr/bin/env bash
# Runs independent training, then fixed minimum-change Procrustes reevaluation.
# Requires legal local DLPFC .h5ad inputs and working GPU software stack.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OUT="${OUT:-${PWD}/LGA_FULL9_OUTPUT}"
export DATA_DIR="${DATA_DIR:?Set DATA_DIR to directory containing dlpfc_<section>.h5ad}"
export BASE="${BASE:-${PWD}}"
export BIB_FULL9_ROOT="$OUT"
export BIB_CANONICAL_OUT="${CANONICAL_OUT:-$OUT/canonical}"
export BIB_FIG6_OUT="${FIG6_OUT:-$OUT/fig6}"
mkdir -p "$OUT" "$BIB_CANONICAL_OUT" "$BIB_FIG6_OUT"
# This first stage uses historical SVD as a replay checkpoint.
bash "$REPO_ROOT/scripts/run_full_9pairs_3seeds.sh"
# This stage asserts the historical results are exactly replayed, then applies
# the geometry-only, pre-fixed minimum-change orthogonal completion.
python3 "$REPO_ROOT/src/current/recompute_canonical.py"
python3 "$REPO_ROOT/src/current/verify_canonical_math.py"
python3 "$REPO_ROOT/src/current/run_fig6_unified.py"
printf '\nOutputs: %s, %s, %s\n' "$OUT" "$BIB_CANONICAL_OUT" "$BIB_FIG6_OUT"
