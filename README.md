# LGA: Latent Geometry Adjustment for Spatial Transcriptomics

**LGA** is a post-hoc approach for improving the compatibility of latent representations constructed from different spatial transcriptomics tissue sections. The method combines anatomical or spatially derived correspondence groups with covariance harmonization and orthogonal Procrustes alignment. It operates on section-level embeddings and does not require retraining the representation encoder during geometric adjustment.

This repository accompanies the manuscript:

> Jia, B. and Yang, P. *Latent Coordinate Compatibility in Spatial Transcriptomics: Anatomy-Guided Post-hoc Geometric Harmonization of Section-Specific Representations.*

## Overview

LGA provides two correspondence settings:

- **LGA-Auto** uses spatial-depth groups derived from tissue coordinates without using ground-truth anatomical labels for alignment.
- **LGA-Oracle** uses known anatomical correspondences to characterize the alignment attainable with supplied labels.

The implementation supports conditional covariance whitening, orthogonal alignment, and evaluation of cross-section label transfer. Experiments cover controlled latent rotations, representation-construction differences, and a mouse-brain spatial pseudo-label example.

## Installation

Clone the repository and create the scientific Python environment:

```bash
git clone https://github.com/3157393972-debug/LGA-spatial-transcriptomics.git
cd LGA-spatial-transcriptomics
conda env create -f environment.yml
conda activate lga-bib-paper1
```

The environment records dependency versions used for the research workflow. CUDA/PyTorch dependencies may require platform-specific installation; see [environment details](docs/ENVIRONMENT.md).

## Quick start

Run the minimal example using simulated embeddings and tissue coordinates:

```bash
python examples/basic_usage.py
```

The core API can also be imported from the repository root:

```python
from lga import auto_groups, axis_from_source, harmonize

axis = axis_from_source(source_coordinates)
source_groups, target_groups, _ = auto_groups(
    source_coordinates, target_coordinates, axis, K=7
)
adjusted_source, adjusted_target, diagnostics = harmonize(
    source_embeddings, target_embeddings,
    source_groups, target_groups,
    whitening="conditional", tau=0.1,
    completion="min_change_sensitivity_only"
)
```

Here, the embeddings are NumPy arrays with the same latent dimension, and the coordinates are spot-by-2 arrays. Each correspondence group used for fitting must contain enough spots in both sections (20 by default). The resulting embeddings can be evaluated using the protocol appropriate to the experiment.

## Reproducing numerical results

Check the supplied evaluation tables and regression tests:

```bash
python -m unittest discover -s tests -v
python scripts/verify_aggregates.py
```

To train the independent GraphSAGE/PCA representations and recompute the associated evaluations, prepare the input AnnData files and run:

```bash
DATA_DIR=/path/to/dlpfc OUT="$PWD/outputs/independent" \
  bash scripts/run_independent_and_reanalyse.sh
```

This data-dependent workflow requires the specified raw inputs and scientific dependencies. Offline tests validate stored outputs and mathematical components; they are not substitutes for full retraining. See [reproduction instructions](docs/REPRODUCING_RESULTS.md) and [protocol definitions](docs/PROTOCOL_BOUNDARIES.md).

## Data

The human dorsolateral prefrontal cortex (DLPFC) sections originate from the [spatialLIBD / Maynard et al. dataset](https://research.libd.org/spatialLIBD/). The mouse-brain sections originate from [10x Genomics](https://www.10xgenomics.com/datasets). The repository includes data specifications and SHA-256 fingerprints; the expression matrices must be obtained separately.

Input preparation, identifiers, annotation requirements and sources are described in [Data and provenance](docs/DATA_AND_PERMISSIONS.md).

## Repository structure

| Path | Contents |
| --- | --- |
| `lga/` | Core LGA API |
| `src/current/` | Independent representation training and numerical evaluation |
| `src/historical/` | Implementations corresponding to the original experimental protocols |
| `scripts/` | Experiment entrypoints and validation utilities |
| `examples/` | Minimal API usage |
| `results/` | Selected numerical outputs |
| `tests/` | Regression and consistency tests |
| `docs/` | Data requirements and experimental details |

## Citation

If you use LGA, please cite the associated manuscript:

```bibtex
@unpublished{JiaYangLGA2026,
  author = {Jia, Bojia and Yang, Ping},
  title = {Latent Coordinate Compatibility in Spatial Transcriptomics: Anatomy-Guided Post-hoc Geometric Harmonization of Section-Specific Representations},
  year = {2026},
  note = {Manuscript}
}
```

## License

This project is licensed under the [MIT License](LICENSE). Copyright (c) 2026 Bojia Jia and Ping Yang.
