"""Check all committed Figure 5/6 plot inputs without external scientific dependencies."""
import csv
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "figure_inputs"
VERIFIED = ROOT / "results" / "verified"

class FigureInputTests(unittest.TestCase):
    def test_figure5_hvg_data(self):
        vals = json.loads((DATA / "e1_preprocessing.json").read_text())
        key_to_value = {(str(v["mode"]), int(v["n_hvg"])): float(v["cd"])
                        for v in vals if int(v["n_pca"]) == 50}
        self.assertEqual(len(key_to_value), 6)
        self.assertEqual(set(key_to_value), {
            (mode, k) for mode in ("independent", "shared")
            for k in (1000, 2000, 3000)
        })
        self.assertTrue(all(math.isfinite(x) and x >= 0 for x in key_to_value.values()))

    def test_covariance_encoder_numeric_inputs(self):
        obj = json.loads((DATA / "encoder_sharing.json").read_text())
        rows = obj if isinstance(obj, list) else obj["results"]
        self.assertIsInstance(rows, list)
        self.assertEqual(len(rows), 12)
        seen = set()
        for r in rows:
            key = (r["donor"], r["pair"])
            self.assertNotIn(key, seen)
            seen.add(key)
            self.assertIn(r["donor"], {"Donor1", "Donor2", "Donor3"})
            for side in ("shared", "independent"):
                x = float(r[side]["CD"])
                self.assertTrue(math.isfinite(x) and x >= 0)

    def test_figure6_pair_means_recomputed_from_rotations(self):
        with (VERIFIED / "FIG6_108_ROTATIONS.csv").open(newline="") as f:
            rotated = list(csv.DictReader(f))
        with (VERIFIED / "FIG6_18_PAIR_MEANS.csv").open(newline="") as f:
            pairs = list(csv.DictReader(f))
        self.assertEqual(len(pairs), 18)
        for row in pairs:
            subset = [r for r in rotated if (r["model"], r["donor"], r["pair"])
                      == (row["model"], row["donor"], row["pair"])]
            self.assertEqual(len(subset), int(row["n_rotated"]))
            for k in ("baseline_LT", "rotated_LT", "auto_LT", "delta_LT"):
                observed = sum(float(r[k]) for r in subset) / len(subset)
                self.assertAlmostEqual(observed, float(row[k]), places=10)

if __name__ == "__main__":
    unittest.main()
