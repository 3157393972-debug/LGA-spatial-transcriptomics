# Evaluation protocols

The manuscript evaluates latent-coordinate compatibility under several settings. Each setting uses the correspondence and scoring rules specified for that experiment.

| Analysis | Representation and correspondence | Evaluation |
| --- | --- | --- |
| Controlled perturbation, Figures 2–4 | Paired training followed by orthogonal perturbations; numerical SVD Procrustes completion | Anatomical label transfer across DLPFC sections |
| Representation-construction mismatch, Figure 5 | Shared or section-specific feature selection and encoder parameters | Whole-section covariance distances and correspondence-informed LGA-Oracle |
| Independent-section analyses, Tables S7–S8 | Separately fitted GraphSAGE and PCA representations; minimum-change Procrustes completion | Pair-level and donor-level label transfer and balanced accuracy |
| Main Figure 6 | PCA, random projection, GCN-AE and VAE under controlled perturbation | Source-to-target nearest-neighbor label transfer |
| Mouse dorsal-cortex analysis, Table S5 | Ordered spatial pseudo-groups for fitting | Full-target transfer of those spatial pseudo-groups |

For paired human LT and the S1B group-count analysis, the target set is defined by anatomical classes present in both sections and is held fixed across evaluated conditions.

LGA-Auto forms correspondences from an ordered spatial axis and a tissue-informed group count. LGA-Oracle forms correspondences from anatomical labels. In the mouse analysis, the pseudo-groups also define the scoring labels. Donors, section pairs, seeds and rotations have different replication roles and are identified with their results.

