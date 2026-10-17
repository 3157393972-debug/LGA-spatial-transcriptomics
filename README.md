# LGA: Latent Geometry Adjustment for Spatial Transcriptomics

**Latent Geometry Adjustment (LGA)** is a post-hoc method for harmonizing latent coordinates across spatial transcriptomics sections. Given section-level embeddings and corresponding anatomical or spatial groups, LGA conditionally applies covariance whitening and estimates an orthogonal Procrustes transformation. The representation encoders do not need to be retrained during this adjustment.

## Method

LGA supports two correspondence settings:

- **LGA-Auto:** constructs ordered groups from tissue spatial coordinates without using spot-level anatomical labels for fitting. It nevertheless requires an anatomically meaningful spatial-depth axis and a grouping resolution appropriate for the tissue.
- **LGA-Oracle:** uses supplied anatomical labels to form correspondences. Results evaluated using those same labels are correspondence-informed diagnostics, not independent label-free validation.

The alignment takes source and target embedding matrices with the same latent dimension, plus one group assignment per spot in each section. The implementation provides conditional ZCA whitening and orthogonal Procrustes alignment. See [protocol boundaries](docs/PROTOCOL_BOUNDARIES.md) for differences between the original and independent-representation experiments.

## Installation

```bash
git clone https://github.com/3157393972-debug/LGA-spatial-transcriptomics.git
cd LGA-spatial-transcriptomics
conda env create -f environment.yml
conda activate lga-bib-paper1
```

The environment records versions from the research setup; it is not a validated, platform-independent lockfile. See [environment and dependencies](docs/ENVIRONMENT.md) for details.

## Quick start

Run the complete, synthetic example:

```bash
python examples/basic_usage.py
```

The example imports `auto_groups`, `axis_from_source`, and `harmonize` from `lga`. It builds spatial-depth groups, aligns simulated embeddings, and reports alignment diagnostics. No biological datasets are needed. The example uses the API's default SVD-based Procrustes completion; experiment-specific minimum-change completion is documented separately. This is an API demonstration, not biological validation.

## Reproducing numerical results

Run the offline consistency checks:

```bash
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
```

These checks test the algorithm and verify selected archived numerical outputs; they do **not** retrain models from raw expression matrices.

For the independent GraphSAGE/PCA experiments, prepare the required DLPFC AnnData files, then run:

```bash
DATA_DIR=/path/to/dlpfc OUT="$PWD/outputs/independent" \
  bash scripts/run_independent_and_reanalyse.sh
```

This workflow requires external data and additional scientific dependencies; a clean-environment, end-to-end rerun has not been established. The original controlled-rotation experiments, independent-representation experiments, main Figure 6 and Supplementary Figure S2 use distinct evaluation protocols. See [reproduction instructions](docs/REPRODUCING_RESULTS.md) before comparing their results.

## Data

The human dorsolateral prefrontal cortex sections are from the [Maynard/spatialLIBD dataset](https://research.libd.org/spatialLIBD/), and the mouse-brain sections are from [10x Genomics](https://www.10xgenomics.com/datasets). Raw or processed expression matrices are not redistributed here. The repository provides input identifiers, a local-input manifest and data-access guidance in [Data and provenance](docs/DATA_AND_PERMISSIONS.md).

## Repository structure

| Path | Contents |
| --- | --- |
| `lga/` | LGA public API |
| `src/current/` | Independent embedding experiments and evaluation |
| `src/historical/` | Historical experimental implementations |
| `scripts/` | Numerical checks and experiment entry points |
| `examples/` | Synthetic API example |
| `tests/` | Unit and regression tests |
| `results/` | Selected archived numerical results |
| `data/` | Input manifest and fingerprints |
| `docs/` | Reproduction, data provenance and protocol guidance |

## Citation

This repository accompanies the following unpublished manuscript:

**Bojia Jia and Ping Yang.** *Latent Coordinate Compatibility in Spatial Transcriptomics: Anatomy-Guided Post-hoc Geometric Harmonization of Section-Specific Representations* (2026).

The bibliographic citation and DOI will be updated when publication details are available.

## License

The source code is distributed under the [MIT License](LICENSE). Copyright (c) 2026 Bojia Jia and Ping Yang.
