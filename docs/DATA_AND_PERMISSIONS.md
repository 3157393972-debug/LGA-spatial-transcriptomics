# Data sources and input files

## Human DLPFC sections

The human experiments use twelve dorsolateral prefrontal cortex sections from the [spatialLIBD study](https://research.libd.org/spatialLIBD/), reported by [Maynard et al. (2021)](https://doi.org/10.1038/s41593-020-00787-0). The section IDs are:

| Donor | Sections |
| --- | --- |
| D1 | 151507, 151508, 151509, 151510 |
| D2 | 151669, 151670, 151671, 151672 |
| D3 | 151673, 151674, 151675, 151676 |

The local analysis files follow the naming pattern `dlpfc_<section>.h5ad`. Their sizes, dimensions and SHA-256 fingerprints are recorded in [the input manifest](../data/paper1_input_manifest.json). The study data and layer annotations can be obtained through [spatialLIBD data access](https://research.libd.org/spatialLIBD/reference/fetch_data.html); the corresponding 10x matrices and image resources are listed in the [spatialLIBD raw-data table](https://github.com/LieberInstitute/spatialLIBD#raw-data).

Anatomical annotations use the spatialLIBD `layer_guess_reordered` field. For D3, the 14,243 barcodes retained across the four local input sections have the same layer annotation as the corresponding official records. The 121 additional official spots have an NA layer annotation (28, 38, 26 and 29 by section). The manuscript reports input-object layer counts in Supplementary Table S9.

## Mouse brain sections

The mouse analyses use 10x Genomics Mouse Brain Serial Section 1, [Sagittal-Anterior](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0) and [Sagittal-Posterior](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-1-0). The local analysis files are `mouse_brain/sagittal_anterior.h5ad` and `mouse_brain/sagittal_posterior.h5ad`; their fingerprints are included in the input manifest. The evaluation selects dorsal-cortex spots and defines spatial pseudo-groups for correspondence and scoring.

## Checking local inputs

Place the study inputs under a local data root and compare their SHA-256 hashes and sizes with the manifest:

```bash
python scripts/verify_input_manifest.py --data-root /path/to/data
```

The manuscript and its data sources should be cited alongside any reuse of these analysis files.
