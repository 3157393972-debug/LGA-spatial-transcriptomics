"""Regression checks for the *revised* implementation; no raw biological data needed."""
import importlib.util
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('independent',ROOT/'src/current/paper1_independent_jsonfixed.py')
independent=importlib.util.module_from_spec(spec)
spec.loader.exec_module(independent)

class TestCore(unittest.TestCase):
    def test_completion_is_orthogonal_and_procrustes_optimal(self):
        rng=np.random.default_rng(42)
        cs=rng.normal(size=(7,32))
        ct=rng.normal(size=(7,32))
        r, rank=independent.canonical_rotation(cs,ct)
        u,_,v=np.linalg.svd(cs.T@ct, full_matrices=True)
        default=u@v
        self.assertLessEqual(np.linalg.norm(cs@r-ct),np.linalg.norm(cs@default-ct)+1e-9)
        np.testing.assert_allclose(r.T@r,np.eye(32),atol=1e-10)
        self.assertLessEqual(np.linalg.norm(r-np.eye(32)),np.linalg.norm(default-np.eye(32))+1e-9)
        self.assertLessEqual(rank,7)

    def test_full_target_includes_unpredicted_class(self):
        source=np.eye(2,dtype=float)
        target=np.array([[1.,0.],[0.,1.],[1.,0.]])
        s=np.array(['L1','L2'])
        t=np.array(['L1','L2','L3'])
        result=independent.full_target_lt(source,target,s,t)
        self.assertEqual(result['n_evaluable'],3)
        self.assertAlmostEqual(result['LT'],2/3)
        self.assertAlmostEqual(result['per_layer_recall']['L3'],0.)

    def test_auto_groups_do_not_use_ground_truth(self):
        src=np.stack([np.arange(50),np.linspace(0,3,50)],axis=1)
        dst=np.stack([np.arange(60),np.linspace(2,6,60)],axis=1)
        gs,gt,edges=independent.auto_groups(src,dst,np.array([1.,0.]),K=7)
        self.assertEqual(len(gs),50)
        self.assertEqual(len(gt),60)
        self.assertEqual(len(edges),6)
        self.assertGreaterEqual(gs.min(),0)
        self.assertLess(gs.max(),7)

if __name__=='__main__':unittest.main()
