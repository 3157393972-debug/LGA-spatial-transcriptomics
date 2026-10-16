# Numerical result reproduction (no plotting code)

This minimal code release deliberately omits figure drawing scripts and presentation assets. The published code retains LGA computations, dataset-dependent experimental execution, and tests of numerical results underlying the manuscript.

Offline checks using verified numerical tables:

```bash
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
```

For latent-level reconstruction from separately available saved embeddings, use `src/current/recompute_canonical.py`, `src/current/verify_canonical_math.py`, and `src/current/run_fig6_unified.py` with `BIB_FULL9_ROOT`. Embedding files are not redistributed. For raw H5AD input, see `scripts/run_independent_and_reanalyse.sh` and `docs/REPRODUCING_RESULTS.md`.

The saved CSV and JSON aggregate inputs remain under `results/` for numerical provenance; they are not rendered to images. This does not establish independent raw-data retraining or reproduction of historical Figures 2–4.
