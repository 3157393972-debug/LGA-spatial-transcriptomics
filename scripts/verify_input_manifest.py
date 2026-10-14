#!/usr/bin/env python3
"""Verify SHA-256 and byte size of the 14 local Paper 1 H5AD inputs.

Uses only the standard library, streams files, and never uploads data.
The manifest hashes are from the author's 2026-10-14 input audit.
"""
import argparse
import hashlib
import json
from pathlib import Path

def sha256_stream(path, block_size=4*1024*1024):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        while chunk := stream.read(block_size):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True, help="Root containing j/ and mouse_brain/")
    parser.add_argument("--manifest", type=Path, default=Path(__file__).resolve().parents[1] / "data" / "paper1_input_manifest.json")
    parser.add_argument("--metadata-only", action="store_true", help="Validate manifest metadata; skip reading raw H5AD inputs")
    args = parser.parse_args()
    obj = json.loads(args.manifest.read_text(encoding="utf-8"))
    rows = obj["files"]
    expected = {str(x) for x in [151507,151508,151509,151510,151669,151670,151671,151672,151673,151674,151675,151676]}
    assert len(rows) == 14, "Expected exactly 14 inputs"
    assert {r["id"] for r in rows} == expected | {"sagittal_anterior","sagittal_posterior"}, "Unexpected input IDs"
    issues = []
    for row in rows:
        rel = Path(row["relative_input_path"])
        if rel.is_absolute() or ".." in rel.parts:
            issues.append(f"Unsafe path in manifest: {rel}")
            continue
        expected_sha = row["sha256"]
        if len(expected_sha) != 64 or any(c not in "0123456789abcdef" for c in expected_sha):
            issues.append(f"Invalid SHA-256 in manifest: {rel}")
            continue
        if args.metadata_only:
            print(f"MANIFEST OK {rel}: {row['n_spots']} spots, {row['n_genes']} genes")
            continue
        path = args.data_root / rel
        if not path.is_file():
            issues.append(f"Missing: {path}")
            continue
        if path.stat().st_size != row["size_bytes"]:
            issues.append(f"Size mismatch: {path}")
            continue
        observed = sha256_stream(path)
        if observed != expected_sha:
            issues.append(f"SHA-256 mismatch: {path}: got {observed}, expected {expected_sha}")
        else:
            print(f"PASS {rel}")
    if issues:
        for problem in issues:
            print("FAIL", problem)
        raise SystemExit(1)
    print(f"Verified {len(rows)} {'manifest records' if args.metadata_only else 'H5AD SHA-256 fingerprints'}.")
if __name__ == "__main__":
    main()
