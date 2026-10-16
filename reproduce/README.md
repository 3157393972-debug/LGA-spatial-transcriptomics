# Minimal figure and result reproduction

The journal-facing reproduction surface keeps only the essential offline checks and Figure 5/6 rendering. Figure 1 is a schematic; its drawing assets are not a runnable science result.

```bash
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
python reproduce/make_figures.py
```

The original analysis script is `src/current/make_figures.py`; this is a compatibility entrypoint. It writes `figures_generated/fig5_a.{pdf,eps}`, `fig5_b.{pdf,eps}` and `fig6.{pdf,eps}`. Input aggregate tables are committed under `results/`. **This reproduces figures from saved outputs, not the training or original Figures 2–4.**

For latent-level reconstruction from archived NPZ, use `src/current/recompute_canonical.py`, `src/current/verify_canonical_math.py` and `src/current/run_fig6_unified.py` with `BIB_FULL9_ROOT`. The NPZs are not redistributed.

All graph generation and other historical paths remain in place for traceability.
