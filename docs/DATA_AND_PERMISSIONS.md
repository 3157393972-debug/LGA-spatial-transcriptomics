# Data access and redistribution review

- Human DLPFC sections (12): 151507–151510, 151669–151672, 151673–151676. The experiment reads local `dlpfc_<section>.h5ad` files. The files are NOT included in this repository candidate. The source dataset, terms of use, annotations, processing and exact file checksums must be documented by the authors before release.
- Mouse anterior/posterior sagittal data: exact preprocessing inputs are `mouse_brain/sagittal_anterior.h5ad` and `mouse_brain/sagittal_posterior.h5ad`; NOT included. No independent anatomic truth labels were demonstrated in the original mouse surrogate-label protocol.
- No raw h5ad, barcodes, private paths, embeddings or trained weights are redistributed in this archive.
- If full 54-record NPZ is later released, verify upstream data licenses and patient-related metadata and record a persistent archive DOI.
- Figure 2–4 historical results used default SVD completion, not the later canonical completion. Results in `results/verified/CANONICAL_ALL_54.csv` use the latter.
