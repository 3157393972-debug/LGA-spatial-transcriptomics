# Reproducing the numerical results

This repository distributes **model, evaluation and numeric verification code only**; figure rendering scripts are not included.

## Offline result checks

```bash
python3 -m compileall -q lga src scripts tests
python3 scripts/verify_aggregates.py
python3 -m unittest discover -s tests -v
```

These verify selected saved numeric summaries, including the 54 independent section/model/seed records and 108 rotated GraphSAGE/PCA records underlying **Supplementary Figure S2**. They do not train models from biological count inputs.

## Supplementary Figure S1B group-count sensitivity

The numerical values are in `results/s1b_group_count_relative.json`. To recalculate them from the nine paired Donor 1 CHNA latent NPZ files saved from the main controlled experiment, run:

```bash
python scripts/recompute_s1b_numeric.py /path/to/latent results/s1b_group_count_relative.json
```

The calculation holds each saved embedding and its ten rotations fixed across group counts. The absolute K=7 mean in this repeat is 0.4460, compared with 0.444 in Table S1 from the original model fits; Supplementary Figure S1B displays the paired within-repeat differences from K=7. The latent NPZ inputs are not redistributed in this repository. This script calculates numerical results only and does not render the figure.

## Independent section-based experiments

The data-dependent workflow requires the 12 local DLPFC `dlpfc_<section>.h5ad` files with spatial coordinates and the expected annotations, plus the training environment:

```bash
DATA_DIR=/path/to/dlpfc OUT="$PWD/outputs/independent" \
 bash scripts/run_independent_and_reanalyse.sh
```

The launcher performs independent training, archived-score verification and fixed minimum-change Procrustes reevaluation. It does **not** generate manuscript figures. It has not been verified in a clean-room GPU rerun.

The separately stored Supplementary Figure S2 rotation results remain in `results/verified/FIG6_108_ROTATIONS.csv`, `FIG6_18_PAIR_MEANS.csv` and `FIG6_DONOR_MEANS.csv` under their **historical filenames**. These files are **not** the final paper's four-representation Figure 6. The numerical CSVs can be verified offline, but this repository does not provide a complete runnable from-input reproduction of those rotations. The filename correspondence is historical, not a scientific equivalence claim.

## Main Figure 6 and original experiments

Main Figure 6 tests PCA, random projection, GCN-AE and VAE in the source-to-target historical protocol. `src/historical/` retains parts of the original cross-method implementations. It does not contain a verified complete runnable implementation for every method in the final four-representation Figure 6; the historical numerical outputs are separately archived. Original Figures 2–4 also use historical implementations. They have not passed full original-input independent reruns; `docs/PROTOCOL_BOUNDARIES.md` explains the differences.

## Mouse experiment

`src/current/mouse_full_target_recheck.py` performs numerical mouse evaluation given external inputs. Its spatial pseudo-groups are used for alignment and scoring, so no independent anatomical ground truth is claimed.

See `docs/DATA_AND_PERMISSIONS.md` and `docs/ENVIRONMENT.md` before attempting raw-data reproduction.
