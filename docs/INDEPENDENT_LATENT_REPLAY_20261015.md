# 2026-10-15 verified independent-latent replay

**Scope:** This is an actual rerun of the analysis on author-provided **saved 32-dimensional independent section-wise embedding NPZs**, not a rerun of GraphSAGE learning on raw `.h5ad` input. The archive `BIB_PAPER1_INDEPENDENT_FULL_9PAIRS.tar.gz` (216 archive members) was staged outside the source repository. No private embeddings or input matrices were committed.

## Executed

A Python 3.13 runtime with available NumPy/SciPy/scikit-learn (not the historically captured Python 3.10.12+Scanpy runtime) was used for the computational replay. `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`.

```bash
export BIB_FULL9_ROOT=/path/to/BIB_PAPER1_INDEPENDENT_FULL_9PAIRS
export BIB_CANONICAL_OUT=/path/to/replay/canonical
export BIB_FIG6_OUT=/path/to/replay/fig6
mkdir -p "$BIB_CANONICAL_OUT" "$BIB_FIG6_OUT"
python src/current/recompute_canonical.py
python src/current/verify_canonical_math.py
python src/current/run_fig6_unified.py
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
```

The scripts now prefer the **versioned repository copy** `src/current/paper1_independent_jsonfixed.py` rather than requiring an external second copy of this Python source. The results still require the original NPZs and 54 pair JSON files under `BIB_FULL9_ROOT`. The sampled source copy in the archive and the staged repository source were byte-identical (SHA-256: `b26b37198b81f689d858c7874f6abea81b4c5c733f47504bceb39df0073a8212`).

## Results from actual latent rerun

| Check | Observed | Saved GitHub results |
|---|---|---|
| Historical default-SVD LT replay in canonical script | All 54 records passed the script's exact score replay assertions | Consistent |
| Canonical 54-row CSV | **0 cell differences** (headers, rows and values) | `results/verified/CANONICAL_ALL_54.csv` |
| Figure 6 108-rotation CSV | **0 cell differences** | `results/verified/FIG6_108_ROTATIONS.csv` |
| Figure 6 six donor-summary rows | **0 cell differences** | `results/verified/FIG6_DONOR_MEANS.csv` |
| Orthogonal completion math | 54/54 records; rank 6 in all; max Frobenius orthogonality error `8.2394e-15`; max absolute objective gap `1.3064e-9`; min-change violation 0 | Passes existing math checks |
| Offline standard-library regression suite | 5 tests passed in approximately 0.002 s | `tests/test_offline_results.py` |

**GraphSAGE:** baseline LT `0.2856850048525696`, canonical LGA-Auto `0.4527659975178706` (27 pair-seed rows, 23 improved).  
**PCA:** baseline `0.43254684193575227`, canonical LGA-Auto `0.4732855797178664` (27 stored rows representing only nine deterministic unique pair outcomes).

**Rotation 108 rows:** GraphSAGE mean rotated LT `0.2484784047117451` -> Auto `0.4512318753627361`; PCA `0.2251719479834889` -> `0.4459039886057816`.

## Current limitations

1. No raw `.h5ad` expression matrices or clean conda installation were used in this replay. Section-wise preprocessing, GraphSAGE training and raw-input reproducibility are **not yet independently rerun**.
2. Figure 5/6 drawing from all Figure 5 JSON inputs in a clean repository checkout is **not yet verified**; Figure 5B data-loader compatibility with both JSON list and `{'results': [...]}` was fixed in source and awaits clean-room renderer execution.
3. Historical Figures 2–4 and mouse full-target rerun have **not** been rebuilt under a clean environment in this round.
4. Author-identified Donor 3 H5AD release provenance and layer annotation origins remain unknown. Source-data rights and dependency/license scans, including all Git history, must be completed before publication.
5. Nothing in this record certifies the repo as public or implies that the raw-data end-to-end reproduction passed. The repository remains private pending author review.

**Version note:** The Python source-file import path was refactored after the local replay, without changing its algorithm or input/output formula. Re-run on the exact final GitHub commit before freezing a public release; this record describes the archived source + equivalent import-path patch checked locally.
