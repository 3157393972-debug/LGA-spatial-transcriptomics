#!/usr/bin/env python3
"""Validate saved spot-level mouse control evidence against JSON diagnostics."""
import argparse, json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--directory',default='/root/autodl-tmp/lga_paper_experiments');a=p.parse_args()
d=Path(a.directory)
data=json.loads((d/'mouse_cortex_sagittal_full_target_checked.json').read_text())
assert len(data['results'])==3

def check(idx, src_labels, target_labels, diagnostic):
    pred=src_labels[idx]
    good=(pred==target_labels)
    common=np.intersect1d(np.unique(target_labels),np.unique(pred))
    keep=np.isin(target_labels,common)&np.isin(pred,common)
    assert len(target_labels)==diagnostic['n_target']
    assert np.sum(keep)==diagnostic['n_masked']
    assert np.sum(~keep)==diagnostic['n_excluded']
    assert np.sum(good)==diagnostic['n_correct']
    assert np.isclose(np.mean(good),diagnostic['lt_full_target'])
    assert np.isclose(np.mean(keep),diagnostic['coverage'])
    assert np.isclose(np.mean(good[keep]) if keep.any() else 0.,diagnostic['lt_masked'])
    recalls=[float(np.mean(good[target_labels==c])) for c in np.unique(target_labels)]
    assert np.isclose(np.mean(recalls),diagnostic['balanced_accuracy_full_target'])

for r in data['results']:
    seed=r['seed'];z=np.load(d/f'mouse_control_seed{seed}_spot_evidence.npz', allow_pickle=False)
    s=z['source_surrogate_labels']; t=z['target_surrogate_labels']
    check(z['original_nn_source_index'],s,t,r['baseline_diagnostic'])
    check(z['unrot_lga_nn_source_index'],s,t,r['unrot_lga_diagnostic'])
    for idx,di in zip(z['rotated_nn_source_indices'],r['rotated_diagnostics']):check(idx,s,t,di)
    for idx,di in zip(z['rotated_lga_nn_source_indices'],r['rotated_lga_diagnostics']):check(idx,s,t,di)
    print(f'seed={seed} n_target={len(t)} full_before={r["lt_orig_full"]:.4f} full_after={r["lt_unrot_lga_full"]:.4f} coverage_before={r["baseline_diagnostic"]["coverage"]:.4f} coverage_after={r["unrot_lga_diagnostic"]["coverage"]:.4f}')
for key,v in data['aggregate'].items():
    assert np.isclose(v,np.mean([r[key] for r in data['results']])),key
print('PASS: all three seeds, full-target metrics and 10 rotations per seed reverified from spot-level prediction indices.')
