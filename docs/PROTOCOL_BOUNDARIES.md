# Evaluation protocols

The manuscript evaluates latent-coordinate compatibility under several settings. Each setting uses the correspondence and scoring rules specified for that experiment.

| Analysis | Representation and correspondence | Evaluation |
| --- | --- | --- |
| Controlled perturbation, Figures 2–4 | Paired training followed by orthogonal perturbations; numerical SVD Procrustes completion | Anatomical label transfer across DLPFC sections |
| Representation-construction mismatch, Figure 5 | Shared or section-specific feature selection and encoder parameters | Whole-section covariance distances and correspondence-informed LGA-Oracle |
| Independent-section analyses, Tables S8–S9 | Separately fitted GraphSAGE and PCA representations; minimum-change Procrustes completion | Pair-level and donor-level label transfer and balanced accuracy |
| Coordinate-selected mouse sagittal analysis, Table S6 | Ordered spatial pseudo-groups for fitting | Full-target transfer of those spatial pseudo-groups |

For paired human LT and the S1B group-count analysis, the target set is defined by anatomical classes present in both sections and is held fixed across evaluated conditions.

LGA-Auto forms correspondences from an ordered spatial axis and a tissue-informed group count. LGA-Oracle forms correspondences from anatomical labels. In the mouse analysis, each section retains spots above the section-specific 60th percentile of the second spatial coordinate. A 10-neighbor graph of spatial coordinates is partitioned by Louvain clustering (resolution 0.3); communities are ordered by mean projection onto the source-derived convex-hull axis and assigned to at most seven pseudo-groups. Those groups define both correspondence and scoring labels. The coordinate selection is not an anatomical boundary. Donors, section pairs, seeds and rotations have different replication roles and are identified with their results.

## Implementation-specific preprocessing and anchors

The archived paired DLPFC controlled-perturbation workflow (`src/historical/unified_experiment_historical.py` and `src/historical/run_rotation_historical.py`) and the mouse full-target workflow (`src/current/mouse_full_target_recheck.py`, which imports the historical implementation) call `lga_harmonize`, which **always** centers and ZCA-whitens both section embeddings before Procrustes. Those historical calls do not apply a covariance-distance threshold. By contrast, the independently fitted GraphSAGE/PCA workflow (`src/current/paper1_independent_jsonfixed.py`) uses `harmonize(..., whitening='conditional', tau=0.1)` for its conditional LGA-Auto result and separately evaluates Procrustes-only and always-ZCA variants. These are distinct experimental implementations, not interchangeable settings.

The upstream paired DLPFC training configuration uses 10-nearest-neighbor mutual-nearest-neighbor (MNN) anchors (`mnn_k=10`), whereas the mouse full-target script calls `generate_anchors_mnn(..., k_nn=20)`. Both use a 10-neighbor **spatial graph**; spatial-graph neighbors and cross-section MNN anchor neighbors are different parameters. The MNN anchors belong to upstream representation training, not the post-hoc LGA group-centroid fitting. Nearest-neighbor label transfer is a downstream evaluation, including baseline and post-alignment comparisons.
