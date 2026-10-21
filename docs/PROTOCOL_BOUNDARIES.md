# Evaluation protocols

The manuscript evaluates latent-coordinate compatibility under several settings. Each setting uses the correspondence and scoring rules specified for that experiment.

| Analysis | Representation and correspondence | Evaluation |
| --- | --- | --- |
| Controlled perturbation, Figures 2–4 | Paired training followed by orthogonal perturbations; numerical SVD Procrustes completion | Anatomical label transfer across DLPFC sections |
| Representation-construction mismatch, Figure 5 | Shared or section-specific feature selection and encoder parameters | Whole-section covariance distances and correspondence-informed LGA-Oracle |
| Independent-section analyses, Tables S7–S8 | Separately fitted GraphSAGE and PCA representations; minimum-change Procrustes completion | Pair-level and donor-level label transfer and balanced accuracy |
| Main Figure 6 | PCA, random projection, GCN-AE and VAE under controlled perturbation | Source-to-target nearest-neighbor label transfer |
| Coordinate-selected mouse sagittal analysis, Table S5 | Ordered spatial pseudo-groups for fitting | Full-target transfer of those spatial pseudo-groups |

For paired human LT and the S1B group-count analysis, the target set is defined by anatomical classes present in both sections and is held fixed across evaluated conditions.

LGA-Auto forms correspondences from an ordered spatial axis and a tissue-informed group count. LGA-Oracle forms correspondences from anatomical labels. In the mouse analysis, each section retains spots above the section-specific 60th percentile of the second spatial coordinate. A 10-neighbor graph of spatial coordinates is partitioned by Louvain clustering (resolution 0.3); communities are ordered by mean projection onto the source-derived convex-hull axis and assigned to at most seven pseudo-groups. Those groups define both correspondence and scoring labels. The coordinate selection is not an anatomical boundary. Donors, section pairs, seeds and rotations have different replication roles and are identified with their results.

