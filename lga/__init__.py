"""Public-facing Paper 1 post-hoc LGA algorithm entrypoint.

Core implementations remain in src/current for exact historical provenance.
No gene-expression data, sample files or plotting packages are imported here.
"""
from .core import auto_groups, axis_from_source, canonical_rotation, covariance_distance, harmonize, zca

__all__ = ["auto_groups", "axis_from_source", "canonical_rotation", "covariance_distance", "harmonize", "zca"]
