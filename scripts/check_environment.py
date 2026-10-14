#!/usr/bin/env python3
"""Report installed versions against the observed AutoDL environment.

This is an environment diagnostic, not proof of end-to-end reproducibility.
"""
import argparse
import importlib.metadata as metadata
import platform
import sys

EXPECTED = {
    "numpy": "2.2.6",
    "scipy": "1.15.3",
    "pandas": "2.3.3",
    "scikit-learn": "1.7.2",
    "scanpy": "1.11.5",
    "anndata": "0.11.4",
    "scikit-misc": "0.5.2",
    "torch": "2.5.1",
    "torch-geometric": "2.8.0.post1",
    "igraph": "1.0.0",
    "louvain": "0.8.2",
    "leidenalg": "0.12.0",
    "umap-learn": "0.5.12",
    "pynndescent": "0.6.0",
    "numba": "0.66.0",
    "h5py": "3.16.0",
    "zarr": "2.18.3",
    "matplotlib": "3.10.9",
    "seaborn": "0.13.2",
    "pillow": "12.3.0",
    "natsort": "8.4.0",
    "tqdm": "4.70.0",
    "joblib": "1.5.3",
    "networkx": "3.4.2",
    "pyyaml": "6.0.3",
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strict', action='store_true', help='Return failure on missing or mismatched dependencies')
    args = parser.parse_args()
    problems = []
    py = platform.python_version()
    if py != '3.10.12':
        problems.append(f'python: observed {py}, reference 3.10.12')
    for package, expected in EXPECTED.items():
        try:
            actual = metadata.version(package)
        except metadata.PackageNotFoundError:
            problems.append(f'{package}: missing (expected {expected})')
            continue
        # Torch wheel local build tags may differ between CPU and CUDA installers.
        if actual.split('+')[0] != expected:
            problems.append(f'{package}: observed {actual}, expected {expected}')
    try:
        import torch
        print(f'torch={torch.__version__}; cuda_available={torch.cuda.is_available()}; device_count={torch.cuda.device_count()}')
    except Exception as exc:
        problems.append(f'torch import failed: {exc}')
    print(f'Python: {py}')
    print(f'Dependencies checked: {len(EXPECTED)}; discrepancies: {len(problems)}')
    for problem in problems:
        print('WARNING:', problem)
    if problems and args.strict:
        sys.exit(1)

if __name__ == '__main__':
    main()
