#!/usr/bin/env bash
set -euo pipefail
BASE="${BASE:-/root/autodl-tmp}"
RUN="$BASE/run_rotation.py"
DATA="$BASE/mouse_brain"
PATCH="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/src/historical/run_mouse_cortex_with_unrotated_control.py"
EXPECTED_RUN='8626503becf68771255ef966b532b97eafe60dcba538d836c2281d97528cf5ff'
for file in "$RUN" "$DATA/sagittal_anterior.h5ad" "$DATA/sagittal_posterior.h5ad"; do
  if [ ! -f "$file" ]; then echo "MISSING REQUIRED INPUT: $file" >&2; exit 2; fi
done
actual="$(sha256sum "$RUN" | awk '{print $1}')"
if [ "$actual" != "$EXPECTED_RUN" ]; then
  echo "Historical run_rotation.py SHA256 differs; abort rather than silently compare incompatible histories" >&2
  echo "Expected $EXPECTED_RUN" >&2; echo "Observed $actual" >&2; exit 3
fi
mkdir -p "$BASE/lga_paper_experiments" "$BASE/BIB_MOUSE_CONTROL_EVIDENCE"
# Imported run_rotation.py is resolved in BASE without editing author source.
export PYTHONPATH="$BASE${PYTHONPATH:+:$PYTHONPATH}"
python3 "$PATCH" 2>&1 | tee "$BASE/BIB_MOUSE_CONTROL_EVIDENCE/run.log"
RESULT="$BASE/lga_paper_experiments/mouse_cortex_sagittal_with_unrotated_control.json"
python3 "$(dirname "${BASH_SOURCE[0]}")/validate_mouse_result.py" "$RESULT"
cp "$RESULT" "$BASE/BIB_MOUSE_CONTROL_EVIDENCE/"
cp "$PATCH" "$BASE/BIB_MOUSE_CONTROL_EVIDENCE/"
cp "$RUN" "$BASE/BIB_MOUSE_CONTROL_EVIDENCE/run_rotation_historical.py"
( cd "$BASE/BIB_MOUSE_CONTROL_EVIDENCE" && sha256sum *.json *.py run.log > SHA256SUMS.txt )
tar -czf "$BASE/BIB_MOUSE_CONTROL_EVIDENCE.tar.gz" -C "$BASE" BIB_MOUSE_CONTROL_EVIDENCE
echo "UPLOAD RESULT: $BASE/BIB_MOUSE_CONTROL_EVIDENCE.tar.gz"
