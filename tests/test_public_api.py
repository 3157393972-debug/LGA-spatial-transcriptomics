"""Regression of the minimal import surface, independent of raw H5AD files."""
import unittest
class PublicAPI(unittest.TestCase):
    def test_core_exports(self):
        from lga import auto_groups, axis_from_source, canonical_rotation
        from lga import covariance_distance, harmonize, zca
        for fn in (auto_groups, axis_from_source, canonical_rotation,
                   covariance_distance, harmonize, zca):
            self.assertTrue(callable(fn))
if __name__ == "__main__":
    unittest.main()
