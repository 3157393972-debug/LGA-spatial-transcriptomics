"""Minimal LGA-Auto example using synthetic embeddings and tissue coordinates.

This example demonstrates the API; it is not a biological evaluation.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from lga import auto_groups, axis_from_source, harmonize


def main():
    rng = np.random.default_rng(42)
    x, y = np.meshgrid(np.arange(21), np.arange(20), indexing="ij")
    source_coords = np.column_stack((x.ravel(), y.ravel())).astype(float)
    target_coords = source_coords + np.array([2.0, -1.0])

    n, dimension = len(source_coords), 16
    source = rng.normal(size=(n, dimension))
    rotation, _ = np.linalg.qr(rng.normal(size=(dimension, dimension)))
    target = source @ rotation + 0.05 * rng.normal(size=(n, dimension))

    axis = axis_from_source(source_coords)
    source_groups, target_groups, _ = auto_groups(
        source_coords, target_coords, axis, K=7
    )
    aligned_source, aligned_target, info = harmonize(
        source, target, source_groups, target_groups,
        whitening="conditional", tau=0.1,
    )

    assert aligned_source.shape == aligned_target.shape == (n, dimension)
    print(
        f"LGA completed: {n} spots, {dimension} dimensions, "
        f"{len(info['common_groups'])} correspondence groups, "
        f"whitening={'on' if info['whitened'] else 'off'}"
    )


if __name__ == "__main__":
    main()
