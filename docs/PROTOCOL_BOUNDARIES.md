# Scientific protocol boundaries

1. Original Figures 2–4: joint-training historical model with numerical SVD completion; the historical recovery rates must not be compared as if recomputed with revised completion.
2. Independent section experiments: independent HVG and encoder construction; min-change Procrustes, saved as 54 section-pair/model/seed records. PCA values are deterministic across nominal seed repeats.
3. Figure 5A covariance discrepancy is normalized; Figure 5B unnormalized Frobenius distance is a different quantity.
4. **Main Figure 6** uses PCA, random projection, GCN-AE, and VAE under archived source-to-target, no-ZCA settings. Its experimental history differs from independent GraphSAGE/PCA.
5. An additional archived rotation analysis uses independently constructed GraphSAGE/PCA embeddings, target-to-source LT, conditional ZCA and min-change completion. CSVs labeled `FIG6_108_ROTATIONS.csv` follow an earlier figure-numbering convention; they do not support a figure in the revised manuscript or the main four-representation Figure 6.
6. Mouse spatial pseudo-labels are used for both grouping and assessment. Mouse full-target calculations are exploratory and are not independent anatomical ground truth.
7. Natural mismatch results can use historical class-intersection evaluation masking; treat them as descriptive.
8. Donor/pair/seed observations are nested, not independent biological replicates. Full raw-data reproduction and unresolved input provenance should not be claimed complete.
