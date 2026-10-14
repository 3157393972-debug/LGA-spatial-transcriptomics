# Environment provenance and remaining lockfile work

Original independent GPU result JSON records: Python **3.10.12**, NumPy **2.2.6**, PyTorch **2.5.1+cu124**. These are *observed versions*, not a complete environment specification.

The requirements file lists needed packages but is intentionally unpinned; no tested conda lock or container digest is included. Scanpy, scikit-misc, python-louvain and platform-specific Torch/CUDA may require version-sensitive installation.

For offline checks, use the existing environment and run `python3 scripts/verify_aggregates.py` and `python3 -m unittest discover -s tests -v`. For a public validated GPU environment, record exact output of `python -m pip freeze`, GPU/driver and CUDA versions, then verify from a fresh checkout on the legally obtained input data.
