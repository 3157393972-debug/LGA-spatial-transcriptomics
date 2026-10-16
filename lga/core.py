"""Compatibility facade for LGA's post-hoc geometry.

Source-of-truth implementation: src/current/paper1_independent_jsonfixed.py.
This wrapper keeps the trained-encoder and evaluation protocols unchanged.
"""
from src.current.paper1_independent_jsonfixed import (
    auto_groups, axis_from_source, canonical_rotation,
    covariance_distance, harmonize, zca,
)
