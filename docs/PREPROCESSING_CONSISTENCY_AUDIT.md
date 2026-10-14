# LGA Paper 1 — input preprocessing consistency audit

These **read-only** scripts check input compatibility; they do not run LGA/GraphSAGE and do not modify H5ADs or send biological data to the network. Audit the complete 14-file input manifest first. **No input dataset has been read by the scripts on the author's AutoDL instance yet.**

## A. Audit all 14 H5AD expression matrices, gene order and labels (no GPU)

From the main GitHub repository checkout:

```bash
python3 scripts/audit_preprocessing_inputs.py \
    --data-root /root/autodl-tmp \
    --output /root/autodl-tmp/BIB_INPUT_MANIFEST/paper1_preprocessing_qc.json
```

By default this samples dispersed blocks of stored CSR `.X/data` values, checks the complete ordered gene identifiers, spatial coordinates, and all DLPFC label counts. To examine **every nonzero stored value** (more disk I/O, still no GPU), add `--full-matrix`.

The report includes the fraction of X entries that are negative, nonfinite or noninteger and summaries of sampled spot-wise sums. Unexpected noninteger counts are a strong **warning** for use of `scanpy.pp.highly_variable_genes(flavor='seurat_v3')`, which expects unnormalized counts. Conversely, integral data *alone* cannot prove equivalence to the original published 10x matrix.

## B. Compare Donor 3 H5AD X against any locally available official 10x filtered matrix (no GPU)

```bash
python3 scripts/compare_dlpfc_10x_counts.py \
    --data-root /root/autodl-tmp \
    --official-dir /root/autodl-tmp/j/dlpfc_raw \
    --n-spots 96 \
    --output /root/autodl-tmp/BIB_INPUT_MANIFEST/donor3_raw_count_comparison.json
```

This maps genes (if uniquely identifiable) and barcodes between the *processed* 151673–151676 `.h5ad` and original 10x filtered `.h5` files. It checks a reproducible sample of 96 complete spot expression vectors per available slice. `--n-spots 0` compares **all** retained spots. A report of `official_precursor_not_available_locally` (notably 151673 in the reported directory) is a **missing comparison**, not a failed numerical check. Use only official filenames when supplying new source matrices; do not substitute mislabeled third-party cached data.

## Why this audit matters

The first eight DLPFC `.h5ad` inputs came from a local spatialLIBD SCE export with synthetic spot IDs and 231 missing layer labels in total; the last four are externally preprocessed, already-filtered H5ADs containing `X_scVI` and `X_C2L`, with unknown exact publisher. The Paper 1 **independent** GraphSAGE/PCA loader uses `.X` as input and runs `seurat_v3` HVG, total count normalization, log1p, scaling and per-section PCA. It does **not** read the imported `X_scVI` or `X_C2L` embeddings. All nine plotted ordered section pairs are within-donor, not across source categories.

The scripts cannot establish correctness of spatial coordinate alignment with original expert annotations after the synthetic barcode conversion. That requires original SCE-to-H5AD conversion provenance or an exact R-side reconciliation. Likewise, the last four external H5ADs' third-party release identification remains open even if their expression X matches official 10x counts exactly.

**Important**: The first script's samples and the second script's spot subsamples are explicit. Do not claim that all matrix values or all labels have been independently validated unless the corresponding full checks have actually been run.