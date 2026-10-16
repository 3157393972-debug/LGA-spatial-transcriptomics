# Latent Geometry Adjustment (LGA) for spatial transcriptomics

**Private staging repository — not yet a public, complete end-to-end reproduction release.**  This repository holds research code and selected verified aggregate results for Paper 1, a *Briefings in Bioinformatics* revision. It is distinct from the later `PAPER2_IDENT` project.

## Publication-facing structure (review branch)

The author-facing, minimal entrypoints are:

- `lga/`: **post-hoc geometric LGA core** (thin import facade, no new algorithm, wraps `src/current/paper1_independent_jsonfixed.py`).
- `experiments/`: **training and evaluation instructions**; working scientific implementations stay in `src/current/` and `src/historical/` to preserve numerical provenance.
- `reproduce/`: **numerical result checks only**, without drawing code.
- `tests/`: result and import regression checks.
- `results/verified/` and `results/figure_inputs/`: essential non-private numerical tables.

This is a **compatibility-first reorganization**, not a file relocation. The `src/current/` and `scripts/` legacy entrypoints continue to work. The optional `lga` facade is not a claim that all source-independent workflows are available without the original biological data. Figure drawing scripts and Figure 1's editable schematic source are excluded; verified numerical evidence and experimental programs remain.

```bash
python -m unittest discover -s tests -v
python scripts/verify_aggregates.py
```

The GitHub Actions offline job now installs the packages needed for core numerical and import smoke tests. A **passing real GitHub Actions run** must still be observed before release; the repository is private. See `reproduce/README.md` and `docs/REPRODUCING_RESULTS.md`.

## Code layout

- `src/historical/`: archived scripts used for original analyses. They include historical evaluation conventions and numerical SVD completion; they are kept for provenance, **not** presented as fully equivalent to the revised pipeline.
- `src/current/`: independent per-section PCA/GraphSAGE experiments, deterministic minimum-change orthogonal completion and unified Figure 6 experiments.
- `scripts/`: verified offline summaries/figure entrypoints, documented data-dependent GPU launchers. Historical launchers are in `archive/scripts/`.
- `results/verified/`: audited, selected descriptive CSV summaries, not raw biological input or a complete set of figures/embeddings.
- `docs/`: data-access requirements and methodological scope.

## Quick checks (do not require raw DLPFC files)

```bash
python3 -m compileall -q src scripts
python3 scripts/verify_aggregates.py
python3 -m unittest discover -s tests -v
```

These checks validate syntax and selected summaries; **they are not a new end-to-end training reproduction**.

## Regenerate independent DLPFC experiments (requires source data)

Original DLPFC input slices: 151507–151510, 151669–151672, and 151673–151676 (12 slices). Supply locally prepared `dlpfc_<section>.h5ad` input files with the expected annotation and spatial columns after checking the upstream terms and preprocessing.

```bash
BASE=/your/workspace DATA_DIR=/your/dlpfc \
  OUT=/your/workspace/independent_outputs \
  bash scripts/run_full_9pairs_3seeds.sh
```

`environment.yml` pins 25 paper-specific dependencies to the versions captured in the author-supplied AutoDL environment; **it is provisional and has not been installed in a fresh environment**. See `docs/ENVIRONMENT.md` and run `python3 scripts/check_environment.py` to compare the local setup. The previous `requirements-not-locked.txt` remains only as an unpinned package inventory. The historical pipeline may have additional requirements and environment-specific paths. The complete raw-data Figure 2–6 pipeline has not been verified from a clean checkout. Figure 5/6 numerical aggregate checks work offline; figure-rendering code is intentionally omitted.

## Mouse experiment: interpretation and boundaries

The additional unrotated LGA-Auto mouse comparison has been executed and independently checked on saved spot predictions: 3 seeds, one mouse section pair, 1,330 post-QC target spots, including 304 spots historically excluded by a class-intersection mask. See `results/verified/mouse_full_target_summary.csv`.

**The mouse groups are spatially derived pseudo-labels used both in alignment and scoring. These results do not constitute independent anatomical ground-truth validation.** The mouse `.h5ad` inputs are not redistributed.

## Verified replay status (2026-10-15)

A new set of five offline consistency tests is available:

```bash
python -m unittest discover -s tests -v
python scripts/verify_aggregates.py
```

The author-supplied **saved embedding archive** has now been used to rerun all 54 canonical evaluations and 108 Figure 6 rotated evaluations; the generated CSV files match the repository's verified CSVs **cell for cell**. This is a **latent-level analysis replay**, not re-training from raw count inputs. The technical record, commands and remaining limitations are in [docs/INDEPENDENT_LATENT_REPLAY_20261015.md](docs/INDEPENDENT_LATENT_REPLAY_20261015.md).

The saved Figure 5/6 numerical inputs are retained for traceability; pure figure-rendering scripts are excluded from this minimal release. Do not announce independent raw-data reproduction or public availability before the release gates are satisfied.

## Reproducibility notes

- Original Figures 2–4 used historical numerical SVD completion; revised independent analyses use a fixed minimum-change orthogonal completion. These protocols are not numerically interchangeable.
- Conditional ZCA, centroid construction and transfer direction also vary in some archived analyses. Do not pool results across scripts without checking evaluation definitions.
- Independent PCA repeated under three nominal seeds has identical embeddings and must not be treated as three independently fitted datasets.
- Donor/pair/seed observations are nested; technical repeats are not independent biological replicates.
- Figure 5A and 5B use different covariance discrepancy normalizations; cross-panel values are not directly comparable.

For runnable paths and important limits, see `docs/REPRODUCING_RESULTS.md`, `docs/PROTOCOL_BOUNDARIES.md` and `docs/DATA_AND_PERMISSIONS.md`. The single-step analysis wrapper is `scripts/run_independent_and_reanalyse.sh`; it is a **GPU entrypoint not verified in a clean environment**.

## Before public release

This is an **author-review staging copy**. Before changing visibility to Public: (1) confirm copyright, third-party dependencies, and input-data licenses; (2) test a clean-checkout reproduction of all claimed main figures and supplement; (3) freeze the environment; (4) check for sensitive paths, credentials and personal data; (5) approve a source-code license and manuscript disclosure; and (6) verify the final paper/rebuttal repository link and immutable commit.

No DOI or GitHub publication status is claimed here.
