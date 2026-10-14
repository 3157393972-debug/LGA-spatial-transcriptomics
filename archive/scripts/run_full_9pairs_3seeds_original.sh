#!/usr/bin/env bash
set -euo pipefail
BASE="${BASE:-/root/autodl-tmp}"
DATA_DIR="${DATA_DIR:-${BASE}/j}"
OUT="${OUT:-${BASE}/BIB_PAPER1_INDEPENDENT_FULL_9PAIRS}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT/logs"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python "$SCRIPT_DIR/paper1_independent_jsonfixed.py" --help >/dev/null
pairs=(
 "151507 151508" "151508 151509" "151509 151510"
 "151669 151670" "151670 151671" "151671 151672"
 "151673 151674" "151674 151675" "151675 151676"
)
for entry in "${pairs[@]}"; do
 read -r s t <<< "$entry"
 for sid in "$s" "$t"; do
   if [ ! -f "$DATA_DIR/dlpfc_${sid}.h5ad" ]; then
     echo "Missing input: $DATA_DIR/dlpfc_${sid}.h5ad" >&2; exit 2
   fi
 done
 echo "=== Running $s -> $t, 3 seeds ==="
 python "$SCRIPT_DIR/paper1_independent_jsonfixed.py" \
   --data-dir "$DATA_DIR" --output "$OUT" \
   --source "$s" --target "$t" \
   --seeds 42 123 456 --models pca graphsage \
   --epochs 500 --K 7 --tau 0.1 \
   --cache-latents --save-latents \
   2>&1 | tee "$OUT/logs/${s}_to_${t}.log"
done
python "$SCRIPT_DIR/summarize_independent.py" "$OUT" | tee "$OUT/logs/summary.log"
cp "$SCRIPT_DIR/paper1_independent_jsonfixed.py" "$OUT/"
cp "$SCRIPT_DIR/summarize_independent.py" "$OUT/"
(cd "$OUT" && find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt)
ARCHIVE="${OUT%/}.tar.gz"
tar -czf "$ARCHIVE" -C "$(dirname "$OUT")" "$(basename "$OUT")"
sha256sum "$ARCHIVE"
echo "COMPLETE $ARCHIVE"
