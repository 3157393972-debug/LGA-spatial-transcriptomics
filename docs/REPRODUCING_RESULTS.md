# Analysis workflows

The repository provides the LGA implementation, experiment entry points, input identifiers and numerical results accompanying the manuscript. The workflows below use the environment in [Environment](ENVIRONMENT.md) and the study inputs described in [Data sources](DATA_AND_PERMISSIONS.md).

## Numerical summaries and software tests

```bash
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
```

The aggregate check reads the saved numerical results, including the 54 independent-section records used for Supplementary Tables S7–S8.

## Supplementary Figure S1B group-count analysis

The values for the group-count comparison are recorded in `results/s1b_group_count_relative.json`. Given the nine saved Donor 1 paired embeddings, calculate the numerical comparison with:

```bash
python scripts/recompute_s1b_numeric.py /path/to/donor1_latents results/s1b_group_count_relative.json
```

The analysis uses the same stored embeddings and rotations across the tested group counts and reports differences relative to K = 7.

## Independently fitted sections

Prepare the twelve local DLPFC files as `dlpfc_<section>.h5ad`, then run:

```bash
DATA_DIR=/path/to/dlpfc OUT="$PWD/outputs/independent" \
  bash scripts/run_independent_and_reanalyse.sh
```

This workflow fits section-specific GraphSAGE and PCA representations and evaluates their cross-section label transfer. Pair, donor and seed results are provided in `results/`.

## Other analysis entry points

The paired controlled-perturbation entry point reads `dlpfc_<section>.h5ad` input files and writes donor JSON summaries:

```bash
DATA_DIR=/path/to/dlpfc OUT_DIR="$PWD/outputs/paired" \
  python src/historical/unified_experiment_historical.py
```

The historical script imports its bundled `run_rotation_historical.py` implementation. It requires the listed scientific dependencies and local study inputs; numerical differences can arise across training environments. The saved result tables are provided separately under `results/`.
- `src/current/mouse_full_target_recheck.py` evaluates mouse full-target spatial pseudo-group transfer using local inputs.
- The main Figure 6 evaluates PCA, random projection, GCN-AE and VAE representations under its source-to-target transfer protocol.

See [Evaluation protocols](PROTOCOL_BOUNDARIES.md) for the training, correspondence and scoring settings used by each analysis.
