import sys, os, time, json, gc
import numpy as np
import scanpy as sc
import torch
import importlib.util
from pathlib import Path
repo_root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('run_rotation', repo_root/'src/historical/run_rotation_historical.py')
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)
from sklearn.preprocessing import normalize

DATA_DIR = os.environ.get('LGA_MOUSE_DATA_DIR','/root/autodl-tmp/mouse_brain')
OUT_DIR = os.environ.get('LGA_MOUSE_OUT','/root/autodl-tmp/lga_paper_experiments')
os.makedirs(OUT_DIR, exist_ok=True)
seeds = [42, 123, 456]

def filter_nan(adata):
    spatial = adata.obsm['spatial']
    mask = ~np.isnan(spatial).any(axis=1)
    return adata[mask].copy()

def fit_depth_hull_normal(adata):
    from scipy.spatial import ConvexHull
    coords = adata.obsm['spatial']
    hull = ConvexHull(coords)
    center = coords.mean(axis=0)
    vertices = coords[hull.vertices]
    farthest_idx = np.argmax(np.linalg.norm(vertices - center, axis=1))
    axis = vertices[farthest_idx] - center
    axis = axis / np.linalg.norm(axis)
    return axis

def auto_label_louvain(adata, axis, resolution=0.3):
    coords = adata.obsm['spatial']
    depth = coords @ axis
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='spatial')
    sc.tl.louvain(adata, resolution=resolution, key_added='louvain_cluster')
    clusters = adata.obs['louvain_cluster'].unique()
    cluster_depth = {c: depth[adata.obs['louvain_cluster'] == c].mean() for c in clusters}
    sorted_clusters = sorted(cluster_depth, key=cluster_depth.get)
    n_clusters = len(sorted_clusters)
    n_layers = min(7, n_clusters)
    mapping = {}
    for i, c in enumerate(sorted_clusters):
        if i < n_layers:
            bucket = int(i * n_layers / n_clusters)
            mapping[c] = bucket
        else:
            mapping[c] = n_layers - 1
    adata.obs['layer_guess'] = adata.obs['louvain_cluster'].map(mapping)
    adata.obs['layer_guess_reordered'] = adata.obs['layer_guess'].apply(lambda x: f'Layer_{x+1}')
    return adata


def evaluate_lt_both(H_s, H_t, labels_s, labels_t):
    """Evaluate target->source nearest neighbor transfer on full target and historical mask.

    True labels here are spatial surrogate group IDs, NOT anatomical ground truth.
    No model fitting or label selection occurs inside the metric.
    """
    assert H_s.shape[0] == len(labels_s) and H_t.shape[0] == len(labels_t)
    H_s_n = normalize(H_s, norm='l2', axis=1)
    H_t_n = normalize(H_t, norm='l2', axis=1)
    nn_idx = (H_t_n @ H_s_n.T).argmax(axis=1)
    predicted = labels_s[nn_idx]
    correct = (predicted == labels_t)
    common = np.intersect1d(np.unique(labels_t), np.unique(predicted))
    mask = np.isin(labels_t, common) & np.isin(predicted, common)
    masked_acc = float(correct[mask].mean()) if mask.any() else 0.0
    full_acc = float(correct.mean())
    per_class = {str(k): float(correct[labels_t == k].mean()) for k in np.unique(labels_t)}
    ret = {
        'lt_masked': masked_acc,
        'lt_full_target': full_acc,
        'balanced_accuracy_full_target': float(np.mean(list(per_class.values()))),
        'per_class_recall': per_class,
        'n_target': int(len(labels_t)),
        'n_masked': int(mask.sum()),
        'n_excluded': int((~mask).sum()),
        'coverage': float(mask.mean()),
        'n_correct': int(correct.sum()),
    }
    if not np.isclose(masked_acc * ret['coverage'], full_acc, atol=1e-12):
        raise AssertionError('Mask excluded a correctly classified spot; inspect evaluation')
    return ret, nn_idx

print('='*70)
print('MOUSE CORTEX (DORSAL SAGITTAL) LGA-AUTO EXPERIMENT')
print('='*70)

# Load
adata_a = filter_nan(sc.read_h5ad(f'{DATA_DIR}/sagittal_anterior.h5ad'))
adata_p = filter_nan(sc.read_h5ad(f'{DATA_DIR}/sagittal_posterior.h5ad'))

# Crop to dorsal cortex (top 40% of y)
y_thresh_a = np.percentile(adata_a.obsm['spatial'][:,1], 60)
y_thresh_p = np.percentile(adata_p.obsm['spatial'][:,1], 60)
cortex_a = adata_a[adata_a.obsm['spatial'][:,1] > y_thresh_a].copy()
cortex_p = adata_p[adata_p.obsm['spatial'][:,1] > y_thresh_p].copy()
print(f'Anterior cortex: {cortex_a.shape[0]} spots, Posterior cortex: {cortex_p.shape[0]} spots')

# Depth axis from cortex of anterior
axis_a = fit_depth_hull_normal(cortex_a)
print(f'Depth axis: ({axis_a[0]:.4f}, {axis_a[1]:.4f})')

# Axis consistency
from scipy.spatial import ConvexHull
hull_a = ConvexHull(cortex_a.obsm['spatial'])
hull_p = ConvexHull(cortex_p.obsm['spatial'])
center_a = cortex_a.obsm['spatial'].mean(axis=0)
center_p = cortex_p.obsm['spatial'].mean(axis=0)
farthest_a = np.argmax(np.linalg.norm(cortex_a.obsm['spatial'][hull_a.vertices] - center_a, axis=1))
farthest_p = np.argmax(np.linalg.norm(cortex_p.obsm['spatial'][hull_p.vertices] - center_p, axis=1))
axis_p = cortex_p.obsm['spatial'][hull_p.vertices[farthest_p]] - center_p
axis_p = axis_p / np.linalg.norm(axis_p)
cos_theta = np.abs(np.dot(axis_a, axis_p))
print(f'Axis consistency |cos| = {cos_theta:.4f}')

# Auto-label
cortex_a = auto_label_louvain(cortex_a, axis_a)
cortex_p = auto_label_louvain(cortex_p, axis_a)

for name, ad in [('Anterior', cortex_a), ('Posterior', cortex_p)]:
    counts = ad.obs['layer_guess_reordered'].value_counts().to_dict()
    print(f'  {name}: {counts}')

# Preprocess
print('\nPreprocessing...')
pp_a = R.preprocess_st_adata(cortex_a.copy(), n_hvg=2000, n_pca=50)
pp_p = R.preprocess_st_adata(cortex_p.copy(), n_hvg=2000, n_pca=50)

data_a = R.prepare_single_slice(pp_a, k=10)
data_p = R.prepare_single_slice(pp_p, k=10)
anchors = R.generate_anchors_mnn(data_a, data_p, k_nn=20)

labels_a = np.array(pp_a.obs['layer_guess_reordered'])
labels_p = np.array(pp_p.obs['layer_guess_reordered'])

results = []
for seed in seeds:
    t0 = time.time()
    out = R.train_model_a(data_a, data_p, anchors, 7, seed)
    H_s, H_t = out['H_s'], out['H_t']
    D = H_s.shape[1]

    base, nn_base = evaluate_lt_both(H_s, H_t, labels_a, labels_p)
    H_s_lga, H_t_w = R.lga_harmonize(H_s, labels_a, H_t, labels_p)
    unrot, nn_unrot = evaluate_lt_both(H_s_lga, H_t_w, labels_a, labels_p)
    # Compare with original function for the unrotated control as a regression check.
    assert np.isclose(base['lt_masked'], R.compute_lt(H_s, H_t, labels_a, labels_p))
    assert np.isclose(unrot['lt_masked'], R.compute_lt(H_s_lga, H_t_w, labels_a, labels_p))

    rng = np.random.RandomState(seed * 1000)
    rotated, rotated_lga = [], []
    nn_rotations, nn_rotations_lga = [], []
    for ri in range(10):
        R_mat = R.random_orthogonal_matrix(D, rng)
        H_s_rot = H_s @ R_mat
        r, nn_r = evaluate_lt_both(H_s_rot, H_t, labels_a, labels_p)
        H_s_rot_lga, _ = R.lga_harmonize(H_s_rot, labels_a, H_t, labels_p)
        a, nn_a = evaluate_lt_both(H_s_rot_lga, H_t_w, labels_a, labels_p)
        rotated.append(r)
        rotated_lga.append(a)
        nn_rotations.append(nn_r)
        nn_rotations_lga.append(nn_a)

    # Historical fields are kept for exact comparability, never overwritten.
    row = {'seed': seed,
           'lt_orig': base['lt_masked'],
           'lt_unrot_lga': unrot['lt_masked'],
           'lt_rot': float(np.mean([x['lt_masked'] for x in rotated])),
           'lt_lga': float(np.mean([x['lt_masked'] for x in rotated_lga])),
           'lt_orig_full': base['lt_full_target'],
           'lt_unrot_lga_full': unrot['lt_full_target'],
           'lt_rot_full': float(np.mean([x['lt_full_target'] for x in rotated])),
           'lt_lga_full': float(np.mean([x['lt_full_target'] for x in rotated_lga])),
           'baseline_diagnostic': base,
           'unrot_lga_diagnostic': unrot,
           'rotated_diagnostics': rotated,
           'rotated_lga_diagnostics': rotated_lga}
    results.append(row)
    print(f"  Seed {seed} ({time.time()-t0:.0f}s): "
          f"masked={base['lt_masked']:.4f}->{unrot['lt_masked']:.4f} "
          f"full={base['lt_full_target']:.4f}->{unrot['lt_full_target']:.4f} "
          f"coverage={base['coverage']:.3f}/{unrot['coverage']:.3f}", flush=True)

    evidence_file = f'{OUT_DIR}/mouse_control_seed{seed}_spot_evidence.npz'
    np.savez_compressed(evidence_file,
        H_s=H_s, H_t=H_t, H_s_unrot_lga=H_s_lga, H_t_whitened=H_t_w,
        source_barcodes=np.asarray(pp_a.obs_names.astype(str)),
        target_barcodes=np.asarray(pp_p.obs_names.astype(str)),
        source_surrogate_labels=labels_a.astype(str),
        target_surrogate_labels=labels_p.astype(str),
        original_nn_source_index=nn_base,
        unrot_lga_nn_source_index=nn_unrot,
        rotated_nn_source_indices=np.stack(nn_rotations),
        rotated_lga_nn_source_indices=np.stack(nn_rotations_lga))
    del out; gc.collect(); torch.cuda.empty_cache()

fields = ['lt_orig','lt_unrot_lga','lt_rot','lt_lga','lt_orig_full','lt_unrot_lga_full','lt_rot_full','lt_lga_full']
aggregate = {k: float(np.mean([r[k] for r in results])) for k in fields}
print('\nSUMMARY: ' + ' '.join(f'{k}={v:.4f}' for k,v in aggregate.items()), flush=True)

with open(f'{OUT_DIR}/mouse_cortex_sagittal_full_target_checked.json', 'w') as f:
    json.dump({'dataset': 'Mouse_Cortex_Sagittal',
        'annotation_type': 'spatial surrogate labels; NOT anatomical ground truth',
        'protocol': 'historical training and LGA, with both historical-masked and full-target LT; 10 rotations per seed',
        'source_surrogate_class_counts': {str(k): int(v) for k,v in zip(*np.unique(labels_a, return_counts=True))},
        'target_surrogate_class_counts': {str(k): int(v) for k,v in zip(*np.unique(labels_p, return_counts=True))},
        'results': results, 'aggregate': aggregate}, f, indent=2)
print('DONE.', flush=True)
