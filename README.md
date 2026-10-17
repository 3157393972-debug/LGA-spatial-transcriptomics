# Latent Geometry Adjustment (LGA) for spatial transcriptomics

Private author-review code for the *Briefings in Bioinformatics* LGA manuscript. This repository contains **method implementations, numerical evaluations, selected verified results, and tests**. It intentionally does **not** provide figure-rendering or layout scripts.

## Structure

- `lga/` — post-hoc geometry API; implementation lives in `src/current/paper1_independent_jsonfixed.py`.
- `src/current/` — independent GraphSAGE/PCA training, minimum-change Procrustes reevaluation, numerical checks and mouse evaluation.
- `src/historical/` — scripts for original experiments and the archived four-representation investigation; see protocol restrictions in `docs/PROTOCOL_BOUNDARIES.md`.
- `scripts/` — experiment launchers and offline numerical verification.
- `results/verified/` — selected output tables (not original count matrices or a full standalone dataset).
- `results/figure_inputs/` — archived numeric evidence for covariance experiments, **not plotting code**.
- `tests/` — model/numerical and aggregate consistency checks.

## Verify saved numeric results (no source tissue matrices needed)

```bash
python3 -m compileall -q lga src scripts tests
python3 -m unittest discover -s tests -v
python3 scripts/verify_aggregates.py
```

These checks verify archived outputs and portions of the mathematical implementation; they do not retrain models from original count data.

## Run independent GraphSAGE/PCA analysis (requires original data)

Prepare `dlpfc_<section>.h5ad` for 12 DLPFC sections as documented in `docs/REPRODUCING_RESULTS.md`. Then run:

```bash
DATA_DIR=/path/to/dlpfc OUT="$PWD/outputs/independent" bash scripts/run_independent_and_reanalyse.sh
```

This trains independently and re-evaluates saved embeddings using minimum-change completion. The raw-data path has **not** been validated from a fresh checkout.

## Manuscript protocol boundaries

The main paper's **Figure 6** is the archived PCA / random projection / GCN-AE / VAE cross-representation experiment. The separately generated **GraphSAGE/PCA 108-rotation** results belong to **Supplementary Figure S2**. Their scoring directions, data representations, geometry settings, and seed protocols are not interchangeable. Numerical data for both studies may be retained without distributing graphics-generation code.

Mouse spatial pseudo-label evaluation does not supply independent anatomical ground truth. Historical Figures 2–4 have not passed clean-room reproduction from original inputs.

Read `docs/REPRODUCING_RESULTS.md`, `docs/PROTOCOL_BOUNDARIES.md`, `docs/DATA_AND_PERMISSIONS.md`, and `docs/ENVIRONMENT.md` for scope and dependencies.

**Publication status:** private staging repository; no public or DOI-bearing release is claimed. Third-party input licenses, full clean-room reruns, and credential/history review remain open.
