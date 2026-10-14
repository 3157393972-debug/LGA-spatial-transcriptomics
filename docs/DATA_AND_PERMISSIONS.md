# Data sources, local provenance and redistribution policy — Paper 1

**Status (2026-10-14): official upstream dataset identities verified from public documentation; local conversion histories supplied by the authors; local input bytes, exact conversion script revisions and sample-by-sample identity remain to be independently checked.**

## 1. Human DLPFC: 12 slices from the Maynard/spatialLIBD project

- Study: Maynard KR et al., *Nature Neuroscience* 24, 425–436 (2021), DOI: [10.1038/s41593-020-00787-0](https://doi.org/10.1038/s41593-020-00787-0).
- Processed 12-section, spot-level study data: [spatialLIBD \`fetch_data(type = "sce")\` reference](https://research.libd.org/spatialLIBD/reference/fetch_data.html). The related R/Bioconductor package paper is Pardo B et al., *BMC Genomics* 23, 434 (2022), DOI: [10.1186/s12864-022-08601-w](https://doi.org/10.1186/s12864-022-08601-w).
- The study team's [spatialLIBD raw-data table](https://github.com/LieberInstitute/spatialLIBD#raw-data) lists filtered 10x \`.h5\` and image resources for all 12 slide IDs, and the original study describes public raw files through the Globus \`jhpce#HumanPilot10x\` endpoint.
- **Only these Paper 1 IDs are included:** \`151507,151508,151509,151510,151669,151670,151671,151672,151673,151674,151675,151676\`.
- The local input naming convention used by this repository is \`dlpfc_<ID>.h5ad\`, under the author's local data directory. These H5ADs are **derived local research artifacts**, not filenames of the official downloaded resources.

### Locally reported lineage (to be checked against exact files and scripts)

| Paper 1 subset | Reported local precursor | Reported conversion | Confidence |
|---|---|---|---|
| 151507–151510, 151669–151672 (8 files) | Local \`Human_DLPFC_Visium_sce.Rdata\` (~407 MB) | \`extract_donors_v2.R\` → \`rds_to_csv_v3.R\`; produced around 2026-08-09 | Author-reported on-machine artifact chain; not yet independently audited against official SCE content |
| 151673–151676 (4 files) | Local individual 10x filtered matrix files | Conversion to H5AD around 2026-08-01 | Author-reported; filtered precursors for **151674–151676** specifically identified; the immediate precursor for **151673** has not yet been documented |

The two paths should **not** be described as a single uniform conversion or pre-normalized dataset. In particular, the official processed spatialLIBD SCE and 10x \`filtered_feature_bc_matrix.h5\` may differ in prior filtering, matrix representation, gene indexing and availability of layer labels. Ensure each generated H5AD contains raw counts as expected by the analysis and has the correct \`obs['layer_guess_reordered']\` annotation and \`obsm['spatial']\` coordinate array. Track the layer-label derivation for the four raw-matrix-converted sections; 10x count H5 files alone do not contain the spatialLIBD expert layer labels.

The separately stored \`dlpfc_Br2720_ant/mid/post\` and \`dlpfc_Br2743_ant/mid/post\` files are from a **different study/project** and must not be described as Paper 1's 12-slice DLPFC cohort or used as its input here.

## 2. Mouse Visium sagittal anterior/posterior

- Official 10x datasets: [Mouse Brain Serial Section 1 (Sagittal-Anterior)](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0) and [Mouse Brain Serial Section 1 (Sagittal-Posterior)](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-1-0), Space Ranger **1.1.0**, released 2020-06-23. Both official dataset pages state **CC BY 4.0**.
- Official upstream file stems are \`V1_Mouse_Brain_Sagittal_Anterior\` and \`V1_Mouse_Brain_Sagittal_Posterior\` (filtered-feature matrix and spatial archive), **not** necessarily the same as local renamed filenames.
- Author-reported local precursor files: \`j/mouse_brain/visium_anterior_filtered_feature_bc_matrix.h5\`, the matching anterior spatial archive, and the corresponding posterior files. The author reports conversion through \`load_mouse_brain.py\`/\`sc.read_visium()\` on 2026-08-12, followed by local restructuring on 2026-08-13.
- Exact local Paper 1 input filenames: \`mouse_brain/sagittal_anterior.h5ad\` and \`mouse_brain/sagittal_posterior.h5ad\` (~175 MB each as reported). The H5ADs were **generated locally** from public-source input data; their bytes have not yet been cross-verified against the official downloadable matrix/spatial files.
- The mouse benchmark uses **spatially derived pseudo-labels** for alignment and evaluation. It does **not** establish independent ground-truth anatomical transfer accuracy.

## 3. Verify exact local study inputs without redistributing raw data

On the author's AutoDL machine, from any shell:

\`\`\`bash
for id in 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676; do
  sha256sum "/root/autodl-tmp/j/dlpfc_\${id}.h5ad"
done
sha256sum /root/autodl-tmp/mouse_brain/sagittal_anterior.h5ad \
          /root/autodl-tmp/mouse_brain/sagittal_posterior.h5ad
\`\`\`

Record checksums, precise input dimensions, barcode matching, the presence/source of label and spatial-coordinate columns, R/Python preprocessing script hashes and output QC counts in a **separate validated manifest**. A filename or file modification timestamp alone is **not proof** of exact file identity, dataset version or provenance. Do not claim those hashes have already been checked.

## 4. Permissions and citations

- The 10x mouse **official dataset pages** license those datasets under CC BY 4.0. Attribute 10x and label derived/preprocessed files as adapted if redistributing them.
- The spatialLIBD **software** package is licensed Artistic-2.0. That software license is **not by itself** a license granting unrestricted redistribution of every DLPFC human-data derivative; verify original dataset terms and attribution separately.
- This repository redistributes **neither raw human/mouse .h5ad files, raw matrices, private latents nor patient-related metadata**. Public upstream access plus a data provenance record is preferable to repackaging files without permissions.
- Cite **both** the Maynard DLPFC study and the spatialLIBD software/data-access paper where the SCE was used. Mouse data should cite the relevant 10x official sample pages.
- The original historical Figure 2–4 pipeline uses historical default SVD completion; newer independent-encoder results in \`results/verified/CANONICAL_ALL_54.csv\` use fixed minimum-change completion. These must not be conflated.

## 5. Outstanding author verification

1. Identify the immediate upstream data/annotation merge path for DLPFC sample **151673**.
2. Validate that both DLPFC conversion branches lead to the expected count-matrix/annotation conventions; inspect the exact R scripts and any local gene/barcode join logs.
3. Compute and preserve SHA-256 of all **14** Paper 1 H5AD inputs and a manifest of their dimensions and relevant columns.
4. Reproduce the mouse pseudo-label experiments on a recorded software environment; do not claim independently verified anatomy labels.
5. Check dataset attribution and permissions, then obtain author signoff before any public code/data release.
