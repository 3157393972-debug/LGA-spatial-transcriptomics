#!/usr/bin/env python3
"""
Unified experiment: single consistent dataset for all LGA paper sections.
Auto-labeling: convex hull depth axis + equal-width depth binning.
Produces: LT_rot, Mdiag_rot, LGA-Oracle, LGA-Auto (auto eval), LGA-Auto (GT eval),
          LGA-Source, ZCA-only, component-wise decomposition.
All 3 donors, 10 rotations per model, 3 seeds.
"""
import warnings; warnings.filterwarnings('ignore')
import sys, os, gc, json, time
import numpy as np
import torch
import scanpy as sc
from scipy.stats import spearmanr
from scipy.spatial import ConvexHull

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_rotation_historical as R

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
DATA_DIR = os.environ.get('DATA_DIR', '/root/autodl-tmp/j')
OUT_DIR = os.environ.get('OUT_DIR', '/root/autodl-tmp/lga_unified')
os.makedirs(OUT_DIR, exist_ok=True)

LAYER_ORDER = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM']
LAYER_NAMES = ['Layer_1', 'Layer_2', 'Layer_3', 'Layer_4', 'Layer_5', 'Layer_6', 'WM']
LAYER_MAP = {'Layer_1': 'L1', 'Layer_2': 'L2', 'Layer_3': 'L3',
             'Layer_4': 'L4', 'Layer_5': 'L5', 'Layer_6': 'L6', 'WM': 'WM'}
LAYER_ORDER_DICT = {'L1': 0, 'L2': 1, 'L3': 2, 'L4': 3, 'L5': 4, 'L6': 5, 'WM': 6}

# ============================================================
# Auto-labeling: convex hull depth axis + depth binning
# ============================================================
def fit_depth_hull_normal(coords):
    hull = ConvexHull(coords)
    vertices = coords[hull.vertices]
    center = coords.mean(axis=0)
    dists = np.linalg.norm(vertices - center, axis=1)
    far_idx = np.argmax(dists)
    axis = vertices[far_idx] - center
    axis = axis / (np.linalg.norm(axis) + 1e-8)
    intercept = -center @ axis
    return axis, intercept

def depth_to_labels(depth_values, n_layers=7):
    depth_norm = (depth_values - depth_values.min()) / (depth_values.max() - depth_values.min() + 1e-8)
    bin_edges = np.linspace(0, 1, n_layers + 1)
    bin_idx = np.digitize(depth_norm, bin_edges[1:-1])
    return np.array([LAYER_NAMES[i] for i in bin_idx]), bin_idx

def auto_label(coords, axis, intercept):
    depth = coords @ axis + intercept
    return depth_to_labels(depth, n_layers=7)

def filter_nan_labels(adata):
    labels = adata.obs['layer_guess_reordered']
    mask = ~labels.isna()
    if mask.sum() < len(labels):
        adata = adata[mask].copy()
    return adata

# ============================================================
# LGA-Source: source GT labels, nearest-centroid for target
# ============================================================
def lga_source(H_s, gt_s, H_t, gt_t):
    H_s_w = R.zca_whiten(H_s)
    H_t_w = R.zca_whiten(H_t)
    ls_m = np.array([LAYER_MAP.get(l, l) for l in gt_s])
    lt_m = np.array([LAYER_MAP.get(l, l) for l in gt_t])
    present = sorted(set(ls_m) & set(lt_m))
    ordered = [l for l in LAYER_ORDER if l in present]
    C_s = np.array([H_s_w[ls_m == l].mean(axis=0) for l in ordered])
    C_s_norm = C_s / (np.linalg.norm(C_s, axis=1, keepdims=True) + 1e-8)
    H_t_w_norm = H_t_w / (np.linalg.norm(H_t_w, axis=1, keepdims=True) + 1e-8)
    tgt_assign = np.argmax(H_t_w_norm @ C_s_norm.T, axis=1)
    C_t = np.array([H_t_w[tgt_assign == k].mean(axis=0) if (tgt_assign == k).sum() > 0
                    else np.zeros(H_t_w.shape[1]) for k in range(len(ordered))])
    U, _, Vt = np.linalg.svd(C_s.T @ C_t)
    H_s_lga = H_s_w @ (U @ Vt)
    return H_s_lga, H_t_w

# ============================================================
# Run one donor
# ============================================================
def run_donor(donor_name, section_ids, pairs):
    print(f"\n{'='*80}")
    print(f"DONOR: {donor_name} | Sections: {section_ids}")
    print(f"{'='*80}", flush=True)

    adatas_raw = {}
    adatas_no_nan = {}
    for sid in section_ids:
        a = sc.read_h5ad(f"{DATA_DIR}/dlpfc_{sid}.h5ad")
        adatas_raw[sid] = a
        adatas_no_nan[sid] = filter_nan_labels(a.copy())
        print(f"  Loaded {sid}: {a.shape[0]} spots", flush=True)

    # Fit depth axis from first section
    first_sid = section_ids[0]
    first_coords = adatas_raw[first_sid].obsm['spatial']
    axis, intercept = fit_depth_hull_normal(first_coords)
    print(f"  Depth axis (from {first_sid}): ({axis[0]:.4f}, {axis[1]:.4f})", flush=True)

    # Auto-label quality
    print(f"  Auto-label quality (depth binning):")
    autolabel_results = {}
    for sid in section_ids:
        coords = adatas_no_nan[sid].obsm['spatial']
        gt = np.array(adatas_no_nan[sid].obs['layer_guess_reordered'])
        pred, _ = auto_label(coords, axis, intercept)
        depth = coords @ axis + intercept
        gt_m = np.array([LAYER_MAP.get(l, l) for l in gt])
        pred_m = np.array([LAYER_MAP.get(l, l) for l in pred])
        gt_o = np.array([LAYER_ORDER_DICT.get(l, -1) for l in gt_m])
        pred_o = np.array([LAYER_ORDER_DICT.get(l, -1) for l in pred_m])
        mask = (gt_o >= 0) & (pred_o >= 0)
        acc = np.mean(gt_o[mask] == pred_o[mask])
        adj = np.mean(np.abs(gt_o[mask] - pred_o[mask]) <= 1)
        r, _ = spearmanr(depth[mask], gt_o[mask])
        print(f"    {sid}: acc={acc:.3f}  adj={adj:.3f}  |r|={abs(r):.3f}", flush=True)
        autolabel_results[sid] = {'acc': acc, 'adj': adj, 'r': abs(r)}

    all_results = []

    for s_id, t_id in pairs:
        print(f"\n  Pair: {s_id} -> {t_id}", flush=True)

        adata_s = R.preprocess_st_adata(adatas_no_nan[s_id].copy(),
            n_hvg=R.CONFIG['n_hvg'], n_pca=R.CONFIG['n_pca'])
        adata_t = R.preprocess_st_adata(adatas_no_nan[t_id].copy(),
            n_hvg=R.CONFIG['n_hvg'], n_pca=R.CONFIG['n_pca'])

        gt_s = np.array(adata_s.obs['layer_guess_reordered'])
        gt_t = np.array(adata_t.obs['layer_guess_reordered'])

        coords_s = adata_s.obsm['spatial']
        coords_t = adata_t.obsm['spatial']

        # Auto-label via depth binning
        auto_s, _ = auto_label(coords_s, axis, intercept)
        auto_t, _ = auto_label(coords_t, axis, intercept)

        data_s = R.prepare_single_slice(adata_s, k=R.CONFIG['spatial_k'])
        data_t = R.prepare_single_slice(adata_t, k=R.CONFIG['spatial_k'])
        anchor_pairs = R.generate_anchors_mnn(data_s, data_t,
            k_nn=R.CONFIG['mnn_k'], top_pct=R.CONFIG['mnn_top_pct'])

        for seed in R.CONFIG['seeds']:
            t0 = time.time()
            print(f"    Training seed {seed}...", end=" ", flush=True)
            outputs = R.train_model_a(data_s, data_t, anchor_pairs,
                                      R.CONFIG['n_communities'], seed)
            elapsed = time.time() - t0
            print(f"({elapsed:.0f}s)", flush=True)

            H_s = outputs['H_s']
            H_t = outputs['H_t']
            D = H_s.shape[1]

            # Baseline
            probe_orig = R.stratified_probe(H_s, gt_s, seed)
            lt_orig = R.compute_lt(H_s, H_t, gt_s, gt_t)
            mdiag_orig = R.compute_mdiag(H_s, H_t, gt_s, gt_t)

            # LGA-Oracle on original
            H_s_lga_orig, H_t_w_orig = R.lga_harmonize(H_s, gt_s, H_t, gt_t)
            lt_lga_orig = R.compute_lt(H_s_lga_orig, H_t_w_orig, gt_s, gt_t)
            mdiag_lga_orig = R.compute_mdiag(H_s_lga_orig, H_t_w_orig, gt_s, gt_t)

            # Rotation experiment
            rot_results = []
            rng = np.random.RandomState(seed * 1000)

            for rot_idx in range(R.CONFIG['n_rotations']):
                R_mat = R.random_orthogonal_matrix(D, rng)
                H_s_rot = H_s @ R_mat

                probe_rot = R.stratified_probe(H_s_rot, gt_s, seed)
                lt_rot = R.compute_lt(H_s_rot, H_t, gt_s, gt_t)
                mdiag_rot = R.compute_mdiag(H_s_rot, H_t, gt_s, gt_t)

                # LGA-Oracle: GT labels
                H_s_rot_lga_gt, H_t_w_gt = R.lga_harmonize(H_s_rot, gt_s, H_t, gt_t)
                lt_rot_lga_gt = R.compute_lt(H_s_rot_lga_gt, H_t_w_gt, gt_s, gt_t)
                mdiag_rot_lga_gt = R.compute_mdiag(H_s_rot_lga_gt, H_t_w_gt, gt_s, gt_t)

                # LGA-Auto (auto eval)
                H_s_rot_lga_auto, H_t_w_auto = R.lga_harmonize(H_s_rot, auto_s, H_t, auto_t)
                lt_rot_lga_auto = R.compute_lt(H_s_rot_lga_auto, H_t_w_auto, auto_s, auto_t)

                # LGA-Auto (GT eval)
                lt_rot_lga_auto_gt = R.compute_lt(H_s_rot_lga_auto, H_t_w_auto, gt_s, gt_t)
                mdiag_rot_lga_auto = R.compute_mdiag(H_s_rot_lga_auto, H_t_w_auto, gt_s, gt_t)

                # LGA-Source
                H_s_rot_lga_src, H_t_w_src = lga_source(H_s_rot, gt_s, H_t, gt_t)
                lt_rot_lga_src = R.compute_lt(H_s_rot_lga_src, H_t_w_src, gt_s, gt_t)

                # ZCA-only
                H_s_rot_zca = R.zca_whiten(H_s_rot)
                H_t_zca = R.zca_whiten(H_t)
                lt_rot_zca = R.compute_lt(H_s_rot_zca, H_t_zca, gt_s, gt_t)

                rot_results.append({
                    'rot_idx': rot_idx,
                    'probe': round(float(probe_rot), 4),
                    'lt': round(float(lt_rot), 4),
                    'mdiag': round(float(mdiag_rot), 4),
                    'lt_lga_oracle': round(float(lt_rot_lga_gt), 4),
                    'mdiag_lga_oracle': round(float(mdiag_rot_lga_gt), 4),
                    'lt_lga_auto': round(float(lt_rot_lga_auto), 4),
                    'lt_lga_auto_gt': round(float(lt_rot_lga_auto_gt), 4),
                    'mdiag_lga_auto': round(float(mdiag_rot_lga_auto), 4),
                    'lt_lga_source': round(float(lt_rot_lga_src), 4),
                    'lt_zca_only': round(float(lt_rot_zca), 4),
                })

            # Aggregate
            def agg(key):
                vals = [r[key] for r in rot_results]
                return round(float(np.mean(vals)), 4), round(float(np.std(vals)), 4)

            probe_rot_m, probe_rot_s = agg('probe')
            lt_rot_m, lt_rot_s = agg('lt')
            mdiag_rot_m, mdiag_rot_s = agg('mdiag')
            lt_lga_oracle_m, lt_lga_oracle_s = agg('lt_lga_oracle')
            mdiag_lga_oracle_m, mdiag_lga_oracle_s = agg('mdiag_lga_oracle')
            lt_lga_auto_m, lt_lga_auto_s = agg('lt_lga_auto')
            lt_lga_auto_gt_m, lt_lga_auto_gt_s = agg('lt_lga_auto_gt')
            mdiag_lga_auto_m, mdiag_lga_auto_s = agg('mdiag_lga_auto')
            lt_lga_source_m, lt_lga_source_s = agg('lt_lga_source')
            lt_zca_m, lt_zca_s = agg('lt_zca_only')

            # Component-wise decomposition
            delta_zca = lt_zca_m - lt_rot_m
            delta_proc = lt_lga_oracle_m - lt_zca_m
            delta_total = lt_lga_oracle_m - lt_rot_m
            if delta_total > 0.001:
                f_zca = delta_zca / delta_total
                f_proc = delta_proc / delta_total
            else:
                f_zca = 0.0; f_proc = 0.0

            print(f"      Probe: orig={probe_orig:.3f} rot={probe_rot_m:.3f}±{probe_rot_s:.3f}")
            print(f"      LT:    orig={lt_orig:.3f} rot={lt_rot_m:.3f}±{lt_rot_s:.3f} "
                  f"ZCA={lt_zca_m:.3f} Oracle={lt_lga_oracle_m:.3f} "
                  f"Auto_gt={lt_lga_auto_gt_m:.3f} Auto={lt_lga_auto_m:.3f} "
                  f"Source={lt_lga_source_m:.3f}")
            print(f"      Mdiag: orig={mdiag_orig:.3f} rot={mdiag_rot_m:.3f}±{mdiag_rot_s:.3f} "
                  f"Oracle={mdiag_lga_oracle_m:.3f} Auto={mdiag_lga_auto_m:.3f}")
            print(f"      Decomp: ZCA={f_zca*100:.1f}% Procrustes={f_proc*100:.1f}%", flush=True)

            all_results.append({
                'donor': donor_name, 'pair': f"{s_id}->{t_id}", 'seed': seed,
                'probe_orig': round(float(probe_orig), 4),
                'lt_orig': round(float(lt_orig), 4),
                'mdiag_orig': round(float(mdiag_orig), 4),
                'lt_lga_orig': round(float(lt_lga_orig), 4),
                'mdiag_lga_orig': round(float(mdiag_lga_orig), 4),
                'probe_rot': probe_rot_m, 'probe_rot_std': probe_rot_s,
                'lt_rot': lt_rot_m, 'lt_rot_std': lt_rot_s,
                'mdiag_rot': mdiag_rot_m, 'mdiag_rot_std': mdiag_rot_s,
                'lt_zca_only': lt_zca_m, 'lt_zca_std': lt_zca_s,
                'lt_lga_oracle': lt_lga_oracle_m, 'lt_lga_oracle_std': lt_lga_oracle_s,
                'mdiag_lga_oracle': mdiag_lga_oracle_m, 'mdiag_lga_oracle_std': mdiag_lga_oracle_s,
                'lt_lga_auto': lt_lga_auto_m, 'lt_lga_auto_std': lt_lga_auto_s,
                'lt_lga_auto_gt': lt_lga_auto_gt_m, 'lt_lga_auto_gt_std': lt_lga_auto_gt_s,
                'mdiag_lga_auto': mdiag_lga_auto_m, 'mdiag_lga_auto_std': mdiag_lga_auto_s,
                'lt_lga_source': lt_lga_source_m, 'lt_lga_source_std': lt_lga_source_s,
                'f_zca': round(f_zca, 4), 'f_proc': round(f_proc, 4),
                'rotations': rot_results, 'time_s': round(float(elapsed), 1),
            })
            del outputs; gc.collect(); torch.cuda.empty_cache()

    # Save
    out_path = f"{OUT_DIR}/{donor_name}_unified.json"
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"  Saved: {out_path}", flush=True)

    # Summary table
    keys = ['lt_orig', 'lt_rot', 'lt_lga_oracle', 'lt_lga_auto_gt', 'lt_lga_auto', 'lt_lga_source',
            'mdiag_orig', 'mdiag_rot', 'mdiag_lga_oracle', 'mdiag_lga_auto', 'f_zca', 'f_proc']
    print(f"\n  {'Pair':<20} {'Seed':>5} ", end='')
    for k in keys:
        print(f"{k:>12}", end=' ')
    print()
    print('-' * (26 + 14 * len(keys)))
    for r in all_results:
        print(f"  {r['pair']:<20} {r['seed']:>5} ", end='')
        for k in keys:
            print(f"{r[k]:>12.4f}", end=' ')
        print()

    print(f"\n  === AGGREGATE {donor_name} ===")
    for k in keys:
        vals = [r[k] for r in all_results]
        print(f"    {k:<25}: {np.mean(vals):.4f} +/- {np.std(vals):.4f}")

    return all_results, autolabel_results


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    t_total = time.time()
    print(f"DEVICE: {DEVICE}", flush=True)
    print(f"UNIFIED LGA EXPERIMENT — depth binning, all 3 donors", flush=True)

    d1_results, d1_al = run_donor('donor1',
        ['151507', '151508', '151509', '151510'],
        [('151507', '151508'), ('151508', '151509'), ('151509', '151510')])

    d2_results, d2_al = run_donor('donor2',
        ['151669', '151670', '151671', '151672'],
        [('151669', '151670'), ('151670', '151671'), ('151671', '151672')])

    d3_results, d3_al = run_donor('donor3',
        ['151673', '151674', '151675', '151676'],
        [('151673', '151674'), ('151674', '151675'), ('151675', '151676')])

    # Cross-donor summary
    print(f"\n{'='*80}")
    print("CROSS-DONOR SUMMARY")
    print(f"{'='*80}")
    print(f"{'Donor':<10} {'LT_orig':>8} {'LT_rot':>8} {'LT_ZCA':>8} "
          f"{'LT_Oracle':>10} {'LT_Auto_gt':>11} {'LT_Auto':>9} {'LT_Source':>10} "
          f"{'Mdiag_rot':>10} {'Mdiag_Ora':>10} {'F_ZCA':>7} {'F_Proc':>7}")
    print('-' * 115)
    for name, results in [('donor1', d1_results), ('donor2', d2_results), ('donor3', d3_results)]:
        m = lambda k: np.mean([r[k] for r in results])
        print(f"{name:<10} {m('lt_orig'):>8.3f} {m('lt_rot'):>8.3f} {m('lt_zca_only'):>8.3f} "
              f"{m('lt_lga_oracle'):>10.3f} {m('lt_lga_auto_gt'):>11.3f} "
              f"{m('lt_lga_auto'):>9.3f} {m('lt_lga_source'):>10.3f} "
              f"{m('mdiag_rot'):>10.3f} {m('mdiag_lga_oracle'):>10.3f} "
              f"{m('f_zca'):>7.3f} {m('f_proc'):>7.3f}")

    # Key findings
    print(f"\n=== KEY FINDINGS ===")
    for name, results in [('donor1', d1_results), ('donor2', d2_results), ('donor3', d3_results)]:
        m = lambda k: np.mean([r[k] for r in results])
        lt_orig = m('lt_orig'); lt_rot = m('lt_rot')
        lt_oracle = m('lt_lga_oracle'); lt_auto_gt = m('lt_lga_auto_gt')
        lt_source = m('lt_lga_source')
        recovery = (lt_auto_gt - lt_rot) / (lt_oracle - lt_rot + 1e-8) * 100
        print(f"  {name}: LT_orig={lt_orig:.3f} LT_rot={lt_rot:.3f} "
              f"Oracle={lt_oracle:.3f} Auto_gt={lt_auto_gt:.3f} "
              f"Source={lt_source:.3f} Recovery={recovery:.0f}%")

    print(f"\nTotal time: {(time.time()-t_total)/60:.1f} min", flush=True)
    print("DONE.", flush=True)
