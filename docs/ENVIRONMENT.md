# Software environment and reproducibility

## Captured AutoDL environment (2026-10-14)

An author-supplied runtime snapshot records **Ubuntu 22.04.1**, **Python 3.10.12**, NumPy **2.2.6**, SciPy **1.15.3**, Scanpy **1.11.5**, AnnData **0.11.4**, scikit-learn **1.7.2**, and PyTorch **2.5.1+cu124**. See `environment.yml` for the pinned first-paper dependency subset.

**Important distinction:** the recorded PyTorch wheel has CUDA 12.4 support, but at the moment this snapshot was captured `torch.cuda.is_available() == False`, the GPU count was zero and `nvidia-smi` returned no device. This proves **no GPU was attached at capture time**; it does *not* independently establish which device was used in every earlier run. A CUDA-enabled wheel does not itself mean training ran on a GPU.

The snapshot reported a 2 GiB container memory limit and 0.5 CPU quota. This is marginal for raw `.h5ad` mouse workflows. A reproduction environment with greater RAM is recommended; exact minimum requirements and end-to-end memory peak have not been independently benchmarked.

## Install from source (provisional, not a solver-tested lock)

```bash
conda env create -f environment.yml
conda activate lga-bib-paper1
python scripts/check_environment.py --strict
python scripts/verify_aggregates.py
python -m unittest discover -s tests -v
```

`environment.yml` pins **25 top-level packages** using versions observed in the supplied `pip freeze`. It is **not** a complete transitive lockfile, and its installation has **not** yet been verified in a fresh container. Installation of PyTorch with CUDA may require a platform-specific wheel index; the YAML pins **2.5.1** instead of the vendor-local `2.5.1+cu124` wheel suffix. If byte-level GPU parity is required, choose the appropriate official PyTorch CUDA 12.4 wheel after recording the GPU/driver and re-run verification. CPU and GPU execution are not guaranteed bitwise-identical.

The full AutoDL `pip freeze` contained approximately 267 distributions, including tooling for other projects. The environment focuses on the packages used by these experiments. This is a minimal **direct-dependency** record, rather than a claim that the whole historical base image has been reproduced.

## Numerical and biological caveats

- The mouse Louvain pseudo-labels may depend sensitively on neighbor-graph implementation, clustering version, spatial-coordinate perturbations, and random state. The author environment note reports such sensitivity; we have not independently reproduced that analysis in a clean environment.
- Mouse pseudo-label agreement is **not** independently annotated anatomical transfer accuracy; changing packages cannot correct that evaluation limitation.
- The original historical Figure 2–4 SVD completion and newer minimum-change completion are different protocols. They remain explicitly separated in source and result provenance.
- `python` was absent as a bare executable in the reported container; activating a conda environment normally supplies `python`. When working directly in the original container, use `python3` if necessary.

## Further reproducibility checks

1. Test `conda env create` and import/smoke tests from a clean Linux installation; record resolved packages (e.g., `conda list --explicit` or an export plus the GPU wheel details).
2. Run the documented end-to-end raw-data workflows on legally obtained inputs and compare Figure 2–6 and mouse evaluations with manuscript outputs.
3. Record data acquisition versions, checksums, figure-specific protocol branches, and machine resources in immutable release notes.
