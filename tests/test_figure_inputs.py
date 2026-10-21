"""Check committed Figure 5 numerical inputs without external scientific dependencies."""
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "figure_inputs"

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


if __name__ == "__main__":
    unittest.main()
