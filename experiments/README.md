# LGA Paper 1 experimental code

The research implementations are maintained under `src/current/` (independent GraphSAGE/PCA training and new results), `src/historical/` (original experiments), and `scripts/` (launchers). These original paths are deliberately **not moved** to preserve archival reproducibility.

From the repository root, run `bash scripts/run_independent_and_reanalyse.sh` after configuring `DATA_DIR`, `OUT` and the dependencies. This requires the verified DLPFC H5AD files and may take considerable time. The full end-to-end training has **not** passed a clean-room rerun.

The new `lga/` facade exposes the post-hoc geometric operators, not a new algorithm implementation. Historical numerical SVD completion and the frozen minimum-change completion are different evaluation protocols. Do not conflate them.

See `docs/REPRODUCING_RESULTS.md`, `docs/INDEPENDENT_LATENT_REPLAY_20261015.md`, and `docs/PUBLIC_RELEASE_REPRODUCTION_GATE_20261015.md`.
