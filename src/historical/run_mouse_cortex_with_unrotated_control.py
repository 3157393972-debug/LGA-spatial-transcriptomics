import sys, os, time, json, gc
import numpy as np
import scanpy as sc
import torch
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_rotation_historical as R

DATA_DIR = '/root/autodl-tmp/mouse_brain'
OUT_DIR = '/root/autodl-tmp/lga_paper_experiments'
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
    
    lt = R.compute_lt(H_s, H_t, labels_a, labels_p)
    H_s_lga, H_t_w = R.lga_harmonize(H_s, labels_a, H_t, labels_p)
    lt_lga = R.compute_lt(H_s_lga, H_t_w, labels_a, labels_p)
    
    rng = np.random.RandomState(seed * 1000)
    lt_rot_vals = []; lt_rot_lga_vals = []
    for ri in range(10):
        R_mat = R.random_orthogonal_matrix(D, rng)
        H_s_rot = H_s @ R_mat
        lt_rot_vals.append(R.compute_lt(H_s_rot, H_t, labels_a, labels_p))
        H_s_rot_lga, _ = R.lga_harmonize(H_s_rot, labels_a, H_t, labels_p)
        lt_rot_lga_vals.append(R.compute_lt(H_s_rot_lga, H_t_w, labels_a, labels_p))
    
    lt_rot = np.mean(lt_rot_vals); lt_rot_lga = np.mean(lt_rot_lga_vals)
    print(f'  Seed {seed} ({time.time()-t0:.0f}s): LT_orig={lt:.3f} LT_unrot_LGA={lt_lga:.3f} LT_rot={lt_rot:.3f} LT_rot_LGA={lt_rot_lga:.3f}', flush=True)
    results.append({'seed': seed, 'lt_orig': float(lt), 'lt_unrot_lga': float(lt_lga), 'lt_rot': float(lt_rot), 'lt_lga': float(lt_rot_lga)})
    del out; gc.collect(); torch.cuda.empty_cache()

avg_orig = np.mean([r['lt_orig'] for r in results])
avg_unrot_lga = np.mean([r['lt_unrot_lga'] for r in results])
avg_rot = np.mean([r['lt_rot'] for r in results])
avg_lga = np.mean([r['lt_lga'] for r in results])
print(f'\nSUMMARY: LT_orig={avg_orig:.3f} LT_unrot_LGA={avg_unrot_lga:.3f} LT_rot={avg_rot:.3f} LT_rot_LGA={avg_lga:.3f}', flush=True)

with open(f'{OUT_DIR}/mouse_cortex_sagittal_with_unrotated_control.json', 'w') as f:
    json.dump({'dataset': 'Mouse_Cortex_Sagittal', 'results': results,
        'aggregate': {'lt_orig': avg_orig, 'lt_unrot_lga': avg_unrot_lga, 'lt_rot': avg_rot, 'lt_lga': avg_lga}}, f, indent=2)
print('DONE.')
