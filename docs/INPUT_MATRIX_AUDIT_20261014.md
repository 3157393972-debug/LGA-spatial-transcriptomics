# Paper 1 input-count and preprocessing audit (2026-10-14)

**Evidence level:** Author-supplied read-only AutoDL execution report `BIB_DATA_AUDIT_REPORT.md` and corresponding `BIB_DATA_AUDIT_SUMMARY.json`, examined for internal consistency on 2026-10-15. The GitHub maintainer did **not** independently download the official matrices or rerun the audit on the author's H5AD input bytes. No private input data or barcode-level tables are redistributed.

## Author-reported data audit

| Test | Reported evidence | Interpretation |
|---|---|---|
| 14 H5AD input matrices | All stored nonzero values integral, nonnegative (fraction integral = 1.0), with variable spot-library sums | Consistent with raw integer UMI counts; no sign of an earlier `normalize_total + log1p` pass |
| Versus official filtered-feature count matrices | 14/14 equal after the report's gene/spot alignment: `nnz(diff)=0`, `max(abs(diff))=0` | Counts on retained spots match the official matrices **according to the executed audit** |
| Eight DLPFC Group A slices | Exact row-fingerprint multiset and matrix comparison; source annotations matched locally exported SCE labels for all nonmissing positions | Gene/count content and labelled rows have strong support; the exported spot IDs themselves are synthetic |
| Four DLPFC Group B slices | 100% of retained rows mapped by count fingerprint; official spots absent: 28 / 38 / 26 / 29 | **Already spot-filtered** official 10x-derived counts; exact third-party H5AD release and layer-label provenance **unknown** |
| Two mouse slices | Barcodes aligned; count matrices reported elementwise equal to 10x sources | Only **input counts** were compared; this does not validate anatomical pseudo-labels or model performance |
| Gene indexing | Gene order unchanged; 24 human and 40 mouse names with Scanpy `-1` uniqueness suffixes | Naming difference, not a count-value difference |
| DLPFC missing layer labels | 231 across Group A; Group B no missing label codes | Must document label mask and denominator in evaluation |
| Input preprocessing | All sections use the same `filter_genes(10) → filter_cells(100) → HVG(seurat_v3, 2000) → normalize_total(1e4) → log1p → scale → PCA(50)` protocol in the historical rotation implementation | No double-normalization indicated; document that later independent representation training may have its own embedding dimensionality |
| Spot QC (DLPFC) | 47,560 input spots → 47,546 after filtering; 14 removed | Aligns with author-provided representation-cache audit |

**Donor naming:** Paper 1 uses Donor 1 = 151507–151510, Donor 2 = 151669–151672, Donor 3 = 151673–151676. Do not use differently named variables from intermediate R scripts as published donor IDs.

## Coordinates and graph geometry: interpret narrowly

The audit records Group A's low-resolution image-pixel coordinates (approximately row/column order) and Group B's full-resolution Space Ranger coordinates (approximately column/row order). The paper's historical `build_spatial_weighted_graph` calls `kneighbors_graph` **separately for each slice** and scales edge lengths with that slice's median nonzero nearest-neighbor distance. This removes a fixed *uniform scale factor* from the formula in exact arithmetic.

The audit's actual numerical tests support **bitwise invariance for an axis swap and a power-of-two rescaling on tested slices**. It also documents numerical differences under translation and substantial kNN edge differences under a 30° rotation because of ties on the hexagonal lattice. Therefore do **not** assert that all rotations, translations, non-binary rescalings, software versions or independently constructed graphs are bitwise invariant.

Separate concern: `fit_depth_hull_normal` chooses an oriented axis geometrically; `depth_to_labels` then min–max normalizes and assigns seven ordered bins. **Min–max normalization does not make the bin labels invariant to reversal of the depth-axis sign:** replacing depth by its negative reverses ordered bin assignments. The audit tested coordinate-system transformation consistency within slices; it does **not** establish that independently chosen pseudo-label bin 1 means the same anatomical layer in different slices. Retain the known caveat about pseudo-label accuracy.

The audit used a revised tissue-image color-saturation mask (`audit5b`) after discarding a degenerate grayscale mask (`audit5`). Tissue overlay and local label autocorrelation are corroborating alignment evidence, **not** independent annotation ground truth for the four Group B sections.

## Integrity review of supplied JSON

Cross-checks of the supplied `BIB_DATA_AUDIT_SUMMARY.json` found:
- All four check sections list the same 14 samples with matching dimensions and QC spot-removal arithmetic.
- 14/14 entries mark `elementwise_exactly_equal=true`, `max_abs_difference=0`, no negative values and integral stored nonzero values.
- The input DLPFC spot total is 47,560; reported QC removes 14; the 231 unlabelled spots are a **different category** and are not the 14 QC-discarded spots.
- The JSON classifies sample **151510** as `square_or_irregular` despite high sixfold angular concentration; hence a blanket “14/14 definitively hex-lattice” label is not fully supported by its own categorical field.

These checks establish **internal consistency of the supplied report and JSON**, not independent verification of upstream binaries or the original audit scripts. The detailed underlying per-file audit JSONs/scripts named in the report were **not** part of the two-file handoff.

## Release and manuscript boundaries

- **Input-count concern substantially resolved:** Donor 3's `adata.X` appears to be raw 10x counts on an already-filtered spot subset, not a log-normalized expression matrix. There is no input-normalization-specific need to rerun the 54 independent-model records.
- **Not yet an end-to-end reproduction:** This audit does not independently reproduce GraphSAGE training, Figures 2–4 or all supplemental figures from a clean checkout.
- The four Group B **expert label origins and exact externally processed object versions** remain unidentified. Report annotation provenance only to the level actually established.
- Do not claim that pseudo-label-based mouse results constitute independent anatomical validation or that depth-axis sign cannot affect matched pseudo-label bins.
- The author report's `recomputation_needed=false` should be understood as “no rerun indicated by this **input counts consistency** check”; it does not close all model, geometry and rebuttal questions.

No raw count matrices, processed H5ADs, cached embeddings, site credentials or intermediate spot-level data are committed to GitHub.
