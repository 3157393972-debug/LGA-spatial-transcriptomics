"""Standard-library-only regression checks on archived LGA Paper 1 CSV results.

These assertions validate saved output consistency, NOT raw-data model reproduction.
"""
import csv
import math
import unittest
from collections import Counter, defaultdict
from pathlib import Path

VERIFIED = Path(__file__).resolve().parents[1] / "results" / "verified"
PAIRS = {
    ("151507", "151508"), ("151508", "151509"), ("151509", "151510"),
    ("151669", "151670"), ("151670", "151671"), ("151671", "151672"),
    ("151673", "151674"), ("151674", "151675"), ("151675", "151676"),
}

def read(name):
    with (VERIFIED / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def mean(rows, key):
    return sum(float(r[key]) for r in rows) / len(rows)


class OfflineResultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.main = read("CANONICAL_ALL_54.csv")

    def test_unique_training_records(self):
        self.assertEqual(len(self.main), 54)
        self.assertEqual(Counter(r["model"] for r in self.main),
                         {"graphsage": 27, "pca": 27})
        keys = {(r["model"], r["source"], r["target"], int(r["seed"]))
                for r in self.main}
        self.assertEqual(len(keys), 54)
        for model in ("graphsage", "pca"):
            for source, target in PAIRS:
                self.assertEqual({seed for m, s, t, seed in keys
                                  if (m, s, t) == (model, source, target)},
                                 {42, 123, 456})

    def test_canonical_values(self):
        expected = {"graphsage": (0.28568500485256954, 0.45276599751787067),
                    "pca": (0.43254684193575227, 0.47328557971786656)}
        for model, (base, auto) in expected.items():
            rows = [r for r in self.main if r["model"] == model]
            self.assertAlmostEqual(mean(rows, "baseline"), base, places=10)
            self.assertAlmostEqual(mean(rows, "Auto_conditional"), auto, places=10)
            for r in rows:
                self.assertGreater(int(r["n_evaluable"]), 0)
                for k in ("baseline", "Auto_conditional", "Auto_procrustes_only",
                          "Oracle_conditional", "baseline_balacc",
                          "Auto_conditional_balacc"):
                    value = float(r[k])
                    self.assertTrue(math.isfinite(value) and 0 <= value <= 1)
                self.assertAlmostEqual(float(r["delta_auto"]),
                                       float(r["Auto_conditional"]) - float(r["baseline"]),
                                       places=11)

    def test_pca_seeds_not_independent(self):
        pairs = defaultdict(list)
        for r in self.main:
            if r["model"] == "pca":
                pairs[(r["source"], r["target"])].append(r)
        self.assertEqual(len(pairs), 9)
        for values in pairs.values():
            for k in ("baseline", "Auto_conditional", "Oracle_conditional"):
                vals = [float(r[k]) for r in values]
                self.assertLess(max(vals) - min(vals), 1e-10)



if __name__ == "__main__":
    unittest.main()
