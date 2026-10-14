# Current release status — 2026-10-14

A 31-file candidate snapshot was uploaded to the **private** repository. Cleanup adds documented analysis entrypoints, unit checks and input summaries. Aggregate/figure regeneration has been checked locally, while raw-data GPU reruns from a fresh checkout remain outstanding.

- Aggregate-only checks and plotting may be run without original input matrices.
- GPU independent-encoder training and historical mouse pipelines **require external raw `.h5ad` files** and matching dependencies; these are not available in this GitHub repository.
- Syntax, aggregate checks, and illustrative unit tests do **not** constitute independent end-to-end reproduction of the original paper.
- Figures 2–4 are supported by historical research scripts, which still require separate data/environment verification.

Do not cite this private staging branch as a final immutable software release or DOI.
