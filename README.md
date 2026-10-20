# LGA: Latent Geometry Adjustment for Spatial Transcriptomics

**Latent Geometry Adjustment (LGA)** harmonizes pre-computed representations across spatial transcriptomics sections. Given section embeddings and corresponding anatomical or ordered spatial groups, LGA conditionally applies covariance whitening and estimates an orthogonal Procrustes transformation.

## Method

LGA supports two correspondence settings:

- **LGA-Auto** constructs ordered groups from a biologically interpretable spatial-depth axis. The group count should reflect the tissue organization.
- **LGA-Oracle** forms group correspondences from supplied anatomical labels. Its results are interpreted as correspondence-informed analyses.

The source and target embeddings must have the same latent dimension, with one group assignment per spot in each section.

## Installation

```bash
git clone https://github.com/3157393972-debug/LGA-spatial-transcriptomics.git
cd LGA-spatial-transcriptomics
conda env create -f environment.yml
conda activate lga-bib-paper1
```

Dependencies and platform setup are described in [Environment](docs/ENVIRONMENT.md).

## Quick start

Run the example on simulated section embeddings:

```bash
python examples/basic_usage.py
```

The example constructs ordered spatial groups, applies LGA and prints alignment diagnostics. The API's SVD-based Procrustes completion and the completion used in specific experiments are described in the [protocol guide](docs/PROTOCOL_BOUNDARIES.md).

## Numerical results and analysis

Check the included numerical summaries and run the software tests:

```bash
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
```

Run the independent GraphSAGE and PCA section analyses with local DLPFC input files:

```bash
DATA_DIR=/path/to/dlpfc OUT="$PWD/outputs/independent" \
  bash scripts/run_independent_and_reanalyse.sh
```

Input identifiers, analysis entry points and the distinct evaluation protocols are described in [Reproducing results](docs/REPRODUCING_RESULTS.md) and [Protocol guide](docs/PROTOCOL_BOUNDARIES.md).

## Data

The human DLPFC sections are available through [spatialLIBD](https://research.libd.org/spatialLIBD/). The mouse-brain sections are available from [10x Genomics](https://www.10xgenomics.com/datasets). Section identifiers and input guidance are listed in [Data and provenance](docs/DATA_AND_PERMISSIONS.md).

## Repository structure

| Path | Contents |
| --- | --- |
| `lga/` | LGA public API |
| `src/current/` | Independent representation experiments and evaluation |
| `src/historical/` | Controlled-perturbation experiment implementations |
| `scripts/` | Numerical checks and analysis entry points |
| `examples/` | Simulated usage example |
| `tests/` | Software and regression tests |
| `results/` | Numerical results |
| `data/` | Input identifiers and fingerprints |
| `docs/` | Analysis protocols, data sources and environment guidance |

## Citation

This repository accompanies the manuscript:

**Bojia Jia and Ping Yang.** *Latent Coordinate Compatibility in Spatial Transcriptomics: Anatomy-Guided Post-hoc Geometric Harmonization of Section-Specific Representations* (2026).

The bibliographic citation and DOI will be updated upon publication.

