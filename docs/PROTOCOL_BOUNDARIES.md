# Scientific provenance, non-equivalence, and negative cases

1. Original Fig. 2–4: historical jointly trained model, default `numpy.linalg.svd` completion; original recovery percentages 81.4/78.3/66.1, not rescored with canonical completion.
2. Independent arm: each DLPFC section performs its own HVG, PCA and/or GraphSAGE, no cross-section MNN or inter-section contrastive loss, minimum-change orthogonal completion; canonical CSV contains 54 pair-seed-model records, PCA replicates identical across three nominal seeds.
3. New Fig. 5A and 5B use covariance discrepancies with DIFFERENT normalizations and thus incomparable absolute units; this replaces an incorrectly described historical PD, not a repaired PD.
4. New Fig. 6 compares two methods, GraphSAGE and PCA, under a common transfer direction and geometric protocol. Earlier four-model claims cannot be inferred from this result.
5. Mouse results use spatial pseudo-labels both for alignment groups and assessment, and DO NOT independently demonstrate anatomy-level prediction.
6. The unrotated+Auto mouse result has been evaluated on all 1,330 post-QC target spots in a separate rerun: full-target LT 0.114286 -> 0.473434. It is a pseudo-label endpoint, not independent anatomical validation; do not combine it with the historical mouse row from a different run.
7. Natural mismatch historical data can retain class-intersection masking and should be treated as exploratory/descriptive.
8. SVD canonical completion is a fixed geometry-only choice, not selected using ground-truth target labels or LT.
9. All reported pair-seed statistics are nested in 9 pairs / 3 donors; do not claim 27 independent biological samples.
