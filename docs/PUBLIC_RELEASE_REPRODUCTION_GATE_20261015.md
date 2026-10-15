## Verified GitHub Actions status (2026-10-15)

Actual PR #1 head `9e9305a6730fcf1d4063121dd5384406a6278045` has a **successful completed** workflow run **38678949265**, `Offline reproducibility checks`. The run's single job `verify-aggregates-and-figures` passed Python 3.10 setup, plotting/core-dependency installation, compileall, archived-aggregate verification, offline `unittest discover`, Figure 5/6 regeneration, output-file checks and artifact upload. [Inspect GitHub Actions run](https://github.com/3157393972-debug/LGA-spatial-transcriptomics/actions/runs/38678949265).

This closes the **real clean-checkout offline/aggregate smoke-test** gap on the PR branch. It does **not** resolve raw-input GraphSAGE training, mouse control, historical Figure 2–4 reproduction, redistribution rights, repository/history security audit, or final license approval. The main branch does not yet include the PR's `lga/`, `reproduce/` and import compatibility tests. Do not mark full scientific end-to-end reproduction complete.

# Paper 1: publication and clean-room reproduction gate (15 Oct 2026)

**Current status: NOT RELEASE-READY. This is a verification plan, not a claim of completed model reruns.** Repository remains private. Do not announce a public reproducibility package until the checks below pass.

## What has actually passed

- `python scripts/verify_aggregates.py` on the author-staged 2026-10-14 source snapshot: GraphSAGE 0.285685 -> 0.452766; PCA 0.432547 -> 0.473286; 108 Figure 6 rotation records. This is **aggregate-file validation**, not a raw-input reproduction.
- `python -m compileall -q src scripts` on that snapshot: Python syntax passes.
- 2026-10-14 author AutoDL audit records all 14 input `adata.X` count matrices as unnormalized integer-valued UMIs and all retained spots as elementwise identical to original 10x count matrices. See [input audit](INPUT_MATRIX_AUDIT_20261014.md) for evidence-level cautions.

## Reproduction gaps found in the staged snapshot

- The staged snapshot lacks `tests/`; `python -m unittest discover -s tests -v` cannot run there. **Do not extrapolate this to current GitHub**, which may include additional files.
- The snapshot `src/current/make_figures.py` expects `src/current/raw_sources/e1_preprocessing.json` and fails in isolation. The **current GitHub version** has been changed to use `results/figure_inputs/e1_preprocessing.json`, `encoder_sharing.json`, `encoder_sharing_cd.json`, and `results/verified/FIG6_18_PAIR_MEANS.csv`. Confirm all four inputs exist in the published tree and the published command succeeds.
- `src/current/run_fig6_unified.py` and `src/current/recompute_canonical.py` still default to an external `BIB_PAPER1_INDEPENDENT_FULL_9PAIRS/paper1_independent_jsonfixed.py` via `BIB_FULL9_ROOT`. The data directory also supplies archived embeddings and expected outputs. These scripts are **not standalone** in a fresh checkout. The exact bundle version, checksum, retrieval instructions and data permissions must be supplied.
- The original Figure 2–4 historical source is still environment-specific. Distinguish historical SVD completion from minimum-change completion used in the revised 54-record experiment.
- Existing `environment.yml` records author versions, but has not been installed/tested from scratch. Inference about reproducibility from pins alone is not sufficient.
- 12 DLPFC inputs + two mouse `.h5ad` objects are excluded intentionally. Raw data access and exact local conversion steps must be documented. Group B labels' original third-party source and provenance remain unresolved.
- Read-only review of the staged snapshot's text files and `.gitignore` found no matches to a small list of conventional secrets; **this is not a comprehensive current-branch + Git-history credential scan**. A credential shared earlier in chat must be rotated, and no credentials should be uploaded.

## Required two-tier independent reproduction

### Level 1: clean checkout, no research inputs
1. Clean clone from an otherwise empty directory.
2. Install/resolve dependencies from a frozen working environment; save `pip freeze`, interpreter version, machine specs.
3. `python -m compileall -q src scripts`
4. `python scripts/verify_aggregates.py`
5. `python -m unittest discover -s tests -v` (only if the published tree has `tests/`; otherwise fix the documentation and add tests).
6. `python src/current/make_figures.py`, compare saved numeric figure inputs and compiled figure exports against the intended revised paper.
7. Archive stdout/stderr, environment fingerprint, code commit and output hashes.

### Level 2: independent training from verified biological input files
1. Verify 14 input file hashes against `data/paper1_input_manifest.json`.
2. On a **fresh directory**, prepare each of 12 DLPFC `.h5ad` according to exact documented provenance; distinguish Group A SCE export from the four Group B externally processed/filtered slices. Verify coordinates, expert labels, raw counts, masking and QC.
3. Independently regenerate per-section HVGs/PCA, separately fitted GraphSAGE embeddings, 9 directed pairs and three GraphSAGE seeds; check each dataset is isolated and no target labels enter Auto/grouping.
4. Re-run and compare **all 54 records** against `results/verified/CANONICAL_ALL_54.csv`; record tolerances in advance (bitwise reproduction may be impossible across hardware/software).
5. Recompute Figure 6's 108 rotated rows, compare the 18 pair means, and re-render Figure 5–6; compare numerical inputs, not only graphical resemblance.
6. Repeat mouse control on the entire 1,330 post-QC target spots; distinguish five source groups versus seven target pseudo groups; no independent ground truth claim.
7. Original Figures 2–4 need their own historical pipeline / environment validation; do not pretend the current standalone train script exactly recreates them.
8. Preserve experiment logs and output checksums. If values differ, reconcile code/data/protocol before claiming independent reproduction.

## Public-release blockers

- Full current-repository file and **Git history** scan for credentials, private paths, unexpected metadata, embedded datasets or patient-related content.
- Rights/license check for third-party code and datasets. **No license has been inferred from the data provider's software license.** Approve a proper source-code license explicitly before attaching one.
- Verify licensing and attribution of any figure assets.
- Merge audited fixes only after review, tag a release with immutable commit and tests, then update manuscript Code/Data Availability *only after* link is accessible.
- User must explicitly switch repository visibility through GitHub settings if no visibility-edit tool is available. A private repository is not publicly reproducible.
- The document is a transparent release-readiness checklist. **It does not authorize changing the study, manuscript or AI-use disclosure.**

### Restrained current claim

“The code repository contains selected archived scripts and verified result aggregates. Checks of the archived aggregates reproduce reported means, but the complete historical and independently trained model pipelines have not yet passed a clean-checkout, raw-input end-to-end reproduction. Public release is pending.”
