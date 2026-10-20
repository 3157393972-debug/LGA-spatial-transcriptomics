# Software environment

The research environment used Python 3.10.12, PyTorch 2.5.1, Scanpy 1.11.5, AnnData 0.11.4, NumPy 2.2.6, SciPy 1.15.3 and scikit-learn 1.7.2 on Ubuntu 22.04. The principal Python dependencies are listed in `environment.yml`.

## Installation

```bash
conda env create -f environment.yml
conda activate lga-bib-paper1
python scripts/check_environment.py --strict
```

For GPU training, install the PyTorch build compatible with the local GPU driver and CUDA runtime.

## Software checks

```bash
python examples/basic_usage.py
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
```

The example exercises the LGA API on simulated embeddings. The other commands check saved numerical summaries and the included tests. See [Analysis workflows](REPRODUCING_RESULTS.md) for data-dependent entry points.

## Evaluation settings

Controlled-perturbation, independently fitted section and mouse analyses use their respective training and scoring procedures. The [evaluation protocol guide](PROTOCOL_BOUNDARIES.md) records these settings, including the Procrustes completion rule and the role of spatial pseudo-groups.
