#!/usr/bin/env python3
"""Compatibility entrypoint for Figure 5 and 6, without duplicating plotting logic."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    runpy.run_path(str(root / "src" / "current" / "make_figures.py"), run_name="__main__")
