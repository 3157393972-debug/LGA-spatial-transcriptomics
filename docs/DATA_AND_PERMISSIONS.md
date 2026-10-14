# Data sources, local provenance and redistribution policy — Paper 1

**Status (2026-10-14): official upstream dataset identities verified from public documentation; local conversion histories supplied by the authors; local input bytes, exact conversion script revisions and sample-by-sample identity remain to be independently checked.**

## 1. Human DLPFC: 12 slices from the Maynard/spatialLIBD project

- Study: Maynard KR et al., *Nature Neuroscience* 24, 425–436 (2021), DOI: [10.1038/s41593-020-00787-0](https://doi.org/10.1038/s41593-020-00787-0).
- Processed 12-section, spot-level study data: [spatialLIBD `fetch_data(type = "sce")` reference](https://research.libd.org/spatialLIBD/reference/fetch_data.html). The related R/Bioconductor package paper is Pardo B et al., *BMC Genomics* 23, 434 (2022), DOI: [10.1186/s12864-022-08601-w](https://doi.org/10.1186/s12864-022-08601-w).
- The study team's [spatialLIBD raw-data table](https://github.com/LieberInstitute/spatialLIBD#raw-data) lists filtered 10x `.h5` and image resources for all 12 slide IDs, and the original study describes public raw files through the Globus `jhpce#HumanPilot10x` endpoint.
- **Only these Paper 1 IDs are included:** `151507,151508,151509,151510,151669,151670,151671,151672,151673,151674,151675,151676`.
- The local input naming convention used by this repository is `dlpfc_<ID>.h5ad`, under the author's local data directory. These H5ADs are **derived local research artifacts**, not filenames of the official downloaded resources.

### Input H5AD fingerprints (author-generated; 2026-10-14)

The 14 exact H5AD inputs now have SHA-256 values, file sizes, spot/gene dimensions and split-provenance flags in [the machine-readable input manifest](../data/paper1_input_manifest.json). The uploaded source is an author-generated read-only AutoDL audit; **the assistant has verified the transcribed entries' structure, not rehashed the 14 raw H5AD bytes independently**.

| Local source group | Slice IDs | H5AD spots | Input provenance |
|---|---|---:|---|
| A (8 files) | 151507–151510, 151669–151672 | 33,317 | Derived in this research workspace from the spatialLIBD SCE, via local R export scripts |
| B (4 files) | 151673–151676 | 14,243 | **Externally preprocessed H5ADs**, original third-party release unidentified; all retained barcodes reported to be a subset of the matching official slices |

Group A has synthetic spot indices and only a `spatial` embedding; the original Visium barcodes were not retained. Group B retains real Visium barcodes suffixed with the sample ID and includes `X_scVI` (50 dimensions), `X_C2L` (75 dimensions), and `uns['CN_order']`. These extra embeddings are not evidence that LGA uses them as model features. The reported Group B barcode audit found 0 nonofficial barcodes, with 28, 38, 26 and 29 official spots missing from 151673, 151674, 151675 and 151676 respectively. This is **a barcode-membership check, not an exact expression-matrix comparison**.

**Correction to earlier provenance inference:** B was *not* generated locally from the original 10x filtered .h5 files. In particular, `151673` has no corresponding precursor in the local `j/dlpfc_raw/` directory. Its exact external processed-data release and the origin of the attached embeddings remain unresolved.

The two paths should **not** be described as a single uniform conversion or pre-normalized dataset. In particular, the official processed spatialLIBD SCE and 10x `filtered_feature_bc_matrix.h5` may differ in prior filtering, matrix representation, gene indexing and availability of layer labels. Ensure each generated H5AD contains raw counts as expected by the analysis and has the correct `obs['layer_guess_reordered']` annotation and `obsm['spatial']` coordinate array. Track the layer-label derivation for the four raw-matrix-converted sections; 10x count H5 files alone do not contain the spatialLIBD expert layer labels.

The separately stored `dlpfc_Br2720_ant/mid/post` and `dlpfc_Br2743_ant/mid/post` files are from a **different study/project** and must not be described as Paper 1's 12-slice DLPFC cohort or used as its input here.

## 2. Mouse Visium sagittal anterior/posterior

- Official 10x datasets: [Mouse Brain Serial Section 1 (Sagittal-Anterior)](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0) and [Mouse Brain Serial Section 1 (Sagittal-Posterior)](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-1-0), Space Ranger **1.1.0**, released 2020-06-23. Both official dataset pages state **CC BY 4.0**.
- Official upstream file stems are `V1_Mouse_Brain_Sagittal_Anterior` and `V1_Mouse_Brain_Sagittal_Posterior` (filtered-feature matrix and spatial archive), **not** necessarily the same as local renamed filenames.
- Author-reported local precursor files: `j/mouse_brain/visium_anterior_filtered_feature_bc_matrix.h5`, the matching anterior spatial archive, and the corresponding posterior files. The author reports conversion through `load_mouse_brain.py`/`sc.read_visium()` on 2026-08-12, followed by local restructuring on 2026-08-13.
- Exact local Paper 1 input filenames: `mouse_brain/sagittal_anterior.h5ad` and `mouse_brain/sagittal_posterior.h5ad` (~175 MB each as reported). The H5ADs were **generated locally** from public-source input data; their bytes have not yet been cross-verified against the official downloadable matrix/spatial files.
- The mouse benchmark uses **spatially derived pseudo-labels** for alignment and evaluation. It does **not** establish independent ground-truth anatomical transfer accuracy.

## 3. Verify exact local study inputs without redistributing raw data

On the author's AutoDL machine, after obtaining the input files at the documented relative paths, check them against the manifest:

```bash
python3 scripts/verify_input_manifest.py --data-root /root/autodl-tmp
```

The verifier streams SHA-256 and checks file size. It **does not** replace inspection of count-matrix semantics, spatial coordinates, cluster/label provenance or the original preprocessing scripts. This work is separate from GPU model reproduction. The author report also records official upstream matrix hashes; these are hashes of **raw .h5 precursors, not processed .h5ad inputs**.

## 4. Permissions and citations

- The 10x mouse **official dataset pages** license those datasets under CC BY 4.0. Attribute 10x and label derived/preprocessed files as adapted if redistributing them.
- The spatialLIBD **software** package is licensed Artistic-2.0. That software license is **not by itself** a license granting unrestricted redistribution of every DLPFC human-data derivative; verify original dataset terms and attribution separately.
- This repository redistributes **neither raw human/mouse .h5ad files, raw matrices, private latents nor patient-related metadata**. Public upstream access plus a data provenance record is preferable to repackaging files without permissions.
- Cite **both** the Maynard DLPFC study and the spatialLIBD software/data-access paper where the SCE was used. Mouse data should cite the relevant 10x official sample pages.
- The original historical Figure 2–4 pipeline uses historical default SVD completion; newer independent-encoder results in `results/verified/CANONICAL_ALL_54.csv` use fixed minimum-change completion. These must not be conflated.

## 5. Outstanding author verification

1. Confirm the exact external processed-data release of DLPFC 151673–151676. Matching official barcodes and gene counts does not prove that expression values and annotations are equivalent.
2. Check the spatialLIBD SCE preprocessing and synthetic barcode reindexing for the first eight sections, and the origin/consistency of layer labels for all twelve.
3. Run the newly added manifest verification against the actual 14 files (without uploading raw files) and retain its output.
4. Test clean-checkout end-to-end model reproduction, including historical and current protocol branches.
5. Secure author approval for source-code and data-use statements prior to publication.
