# Reproducing Paper 1 results

## Offline / no raw datasets (tested separately)

```bash
python3 scripts/verify_aggregates.py
python3 -m unittest discover -s tests -v
```

Saved figure-related numerical inputs remain in `results/figure_inputs/` and `results/verified/`. Pure drawing scripts and generated PDF/EPS figures are excluded from this minimal code release. Figure 5A normalized covariance distance and Figure 5B unnormalized Frobenius covariance distance have incompatible units. Aggregate checks do **not** reproduce original model training.

## Independent GraphSAGE/PCA from raw 12-slice DLPFC dataset (GPU, NOT clean-room tested)

Provide 12 `dlpfc_<slice>.h5ad` with `obsm['spatial']`, count matrix, and annotation column `obs['layer_guess_reordered']` (used for evaluation/Oracle only). Check exact preprocessing and data rights before use. Fix the environment to match manuscript-associated experiments. Then:

```bash
DATA_DIR=/path/to/dlpfc BASE="$PWD" OUT="$PWD/outputs/independent" \
  bash scripts/run_independent_and_reanalyse.sh
```

The launcher (a) trains each section independently and saves NPZ, (b) **replays the historical evaluation** and fails if stored LT scores disagree, (c) applies minimum-change orthogonal completion, (d) generates unified Figure 6 CSV. It needs several GPU hours depending on resources. Results can differ from paper if data preprocessing, model libraries or numerical environments differ. The code currently requires external dependencies not pinned to a validated lockfile.

## Mouse full-target pseudo-label validation (GPU, NOT clean-room tested)

```bash
LGA_MOUSE_DATA_DIR=/path/to/mouse_brain LGA_MOUSE_OUT="$PWD/outputs/mouse" \
  python3 src/current/mouse_full_target_recheck.py
```

This is **pseudo-label agreement**: the spatial groups used for evaluation also enter fitting, so it is not independent histological ground truth. Source/target pseudo-group coverage differs. The complete scored target population is used for full-target accuracy.

## Original Figure 2–4

See `src/historical/` and `docs/PROTOCOL_BOUNDARIES.md`. Older scripts depend on original data layout and do not constitute a one-command reproducible setup. The current minimum-change completion does not retrospectively apply to those historical summaries.

## Data and software rights

See `docs/DATA_AND_PERMISSIONS.md` and repository README. No raw Visium matrices or patient-related embeddings are redistributed here; confirm licenses before public release.
