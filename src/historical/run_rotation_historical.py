#!/usr/bin/env python3
"""
Controlled Rotation Experiment — Core causal evidence for LGA paper.

Design:
  1. Train model with original config → H_s, H_t
  2. Apply random orthogonal rotation R: H_s_rot = H_s @ R
  3. Measure:
     - Probe (within-slice) → should be UNCHANGED (information preserved)
     - LT (cross-slice) → should DROP (coordinate change destroys correspondence)
     - Mdiag → should DROP (coordinate incompatibility)
  4. LGA: H_s_rot → ZCA + Procrustes → H_s_aligned
  5. Measure:
     - LT → should RECOVER
     - Mdiag → should RECOVER

Runs: 3 pairs × 3 seeds × N_rotations → comprehensive evidence
"""
import warnings; warnings.filterwarnings('ignore')
import sys, os, gc, json, time
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import scanpy as sc
from scipy.stats import spearmanr
from sklearn.neighbors import kneighbors_graph
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
DATA_DIR = os.environ.get('DATA_DIR', '/root/autodl-tmp/j')
OUT_DIR = os.environ.get('OUT_DIR', '/root/autodl-tmp/lga_rotation_results')
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# CONFIG (original paper reproduction)
# ============================================================
CONFIG = {
    "n_communities": 7,
    "input_dim": 50, "hidden_dim": 64, "output_dim": 32,
    "temperature": 0.18,
    "n_hvg": 2000, "n_pca": 50,
    "spatial_k": 10, "mnn_k": 10, "mnn_top_pct": 0.05,
    "lambda_anchor": 5.0, "lambda_domain": 0.3, "lambda_community": 0.3,
    "lambda_laplacian": 0.01, "lambda_entropy": 0.03, "lambda_balance": 0.0,
    "lambda_cst_entropy": 0.01, "lambda_mu_align": 0.2,
    "lr": 1e-3, "weight_decay": 1e-4, "n_epochs": 500, "grad_clip": 1.0,
    "seeds": [42, 123, 456],
    "n_rotations": 10,  # number of random rotations per seed
}
LAYER_ORDER = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM']
LAYER_MAP = {'Layer_1': 'L1', 'Layer_2': 'L2', 'Layer_3': 'L3',
             'Layer_4': 'L4', 'Layer_5': 'L5', 'Layer_6': 'L6', 'WM': 'WM'}

# ============================================================
# EXACT ORIGINAL MODEL (model_v6_2.py)
# ============================================================
class GraphSAGELayer(nn.Module):
    def __init__(self, in_dim, out_dim, use_activation=True):
        super().__init__()
        self.linear = nn.Linear(in_dim * 2, out_dim, bias=False)
        self.use_activation = use_activation
    def forward(self, A_norm, H):
        H_cat = torch.cat([H, A_norm @ H], dim=1)
        H_out = self.linear(H_cat)
        return F.leaky_relu(H_out, 0.2) if self.use_activation else H_out

class SpatialEncoder(nn.Module):
    def __init__(self, input_dim=50, hidden_dim=64, output_dim=32):
        super().__init__()
        self.gcn1 = GraphSAGELayer(input_dim, hidden_dim, use_activation=True)
        self.gcn2 = GraphSAGELayer(hidden_dim, output_dim, use_activation=False)
    def forward(self, A_norm, X):
        return self.gcn2(A_norm, self.gcn1(A_norm, X))

class SoftAssignment(nn.Module):
    def __init__(self, temperature=0.18):
        super().__init__()
        self.temperature = temperature
    def forward(self, H, mu):
        return F.softmax(F.normalize(H, dim=1) @ F.normalize(mu, dim=1).T / self.temperature, dim=1)

class CommunityCenter(nn.Module):
    def forward(self, H, R):
        return (R.T @ H) / R.sum(dim=0).clamp(min=1e-12).unsqueeze(1)

class ModularityLoss(nn.Module):
    def forward(self, R, B, total_edge_weight):
        return -(R * (B @ R)).sum() / (total_edge_weight + 1e-12)

class CHNA_ST_Model_Original(nn.Module):
    def __init__(self, input_dim=50, hidden_dim=64, output_dim=32,
                 n_communities_s=7, n_communities_t=None,
                 temperature=0.18, lambda_entropy=0.03, lambda_balance=0.0,
                 lambda_cst_entropy=0.01, lambda_mu_align=0.2):
        super().__init__()
        self.n_communities_s = n_communities_s
        self.n_communities_t = n_communities_t if n_communities_t is not None else n_communities_s
        self.output_dim = output_dim
        self.temperature = temperature
        self.lambda_entropy = lambda_entropy
        self.lambda_balance = lambda_balance
        self.lambda_cst_entropy = lambda_cst_entropy
        self.lambda_mu_align = lambda_mu_align

        self.encoder_s = SpatialEncoder(input_dim, hidden_dim, output_dim)
        self.encoder_t = self.encoder_s
        self.soft_assign = SoftAssignment(temperature=temperature)
        self.center_update = CommunityCenter()
        self.domain_loss_fn = ModularityLoss()

        self.mu_s = nn.Parameter(torch.randn(n_communities_s, output_dim))
        nn.init.xavier_uniform_(self.mu_s)
        self.mu_t = nn.Parameter(torch.randn(self.n_communities_t, output_dim))
        nn.init.xavier_uniform_(self.mu_t)

        K = n_communities_s
        self.C_st_raw = nn.Parameter(torch.eye(K) * 1.0)
        self.C_ts_raw = nn.Parameter(torch.eye(K) * 1.0)

    def get_C_st(self): return F.softmax(self.C_st_raw, dim=1)
    def get_C_ts(self): return F.softmax(self.C_ts_raw, dim=1)

    def forward_dual_slice(self, A_norm_s, X_s, B_s, total_edge_s,
                           A_norm_t, X_t, B_t, total_edge_t, anchor_pairs):
        H_s = self.encoder_s(A_norm_s, X_s)
        H_t = self.encoder_t(A_norm_t, X_t)
        R_s = self.soft_assign(H_s, self.mu_s)
        R_t = self.soft_assign(H_t, self.mu_t)
        mu_s_new = self.center_update(H_s, R_s)
        mu_t_new = self.center_update(H_t, R_t)
        loss_domain_s = self.domain_loss_fn(R_s, B_s, total_edge_s)
        loss_domain_t = self.domain_loss_fn(R_t, B_t, total_edge_t)
        loss_anchor = self._anchor_loss(H_s, H_t, anchor_pairs)
        C_st = self.get_C_st()
        C_ts = self.get_C_ts()
        loss_community = self._community_contrastive_loss(mu_s_new, mu_t_new, C_st, C_ts)
        return {
            "H_s": H_s, "H_t": H_t,
            "R_s": R_s, "R_t": R_t,
            "mu_s": mu_s_new, "mu_t": mu_t_new,
            "loss_domain_s": loss_domain_s, "loss_domain_t": loss_domain_t,
            "loss_anchor": loss_anchor,
            "loss_community": loss_community,
            "loss_entropy_s": self._entropy_loss(R_s),
            "loss_entropy_t": self._entropy_loss(R_t),
            "loss_balance_s": self._balance_loss(R_s),
            "loss_balance_t": self._balance_loss(R_t),
            "loss_cst_entropy": self._cst_entropy_loss(C_st),
            "loss_mu_align": self._mu_alignment_loss(mu_s_new, mu_t_new, C_st),
            "hard_assignment_s": torch.argmax(R_s, dim=1),
            "hard_assignment_t": torch.argmax(R_t, dim=1),
        }

    def _anchor_loss(self, H_s, H_t, anchor_pairs):
        if anchor_pairs.shape[0] == 0: return torch.tensor(0.0, device=H_s.device)
        return (H_s[anchor_pairs[:, 0]] - H_t[anchor_pairs[:, 1]]).norm(p=2, dim=1).mean()

    def _community_contrastive_loss(self, mu_s, mu_t, C_st, C_ts):
        mu_s_n, mu_t_n = F.normalize(mu_s, dim=1), F.normalize(mu_t, dim=1)
        cos = mu_s_n @ mu_t_n.T
        exp_cos = torch.exp(cos)
        denom = exp_cos.sum(dim=1) + 1e-8
        loss_st = -torch.log((C_st * exp_cos).sum(dim=1) / denom).mean()
        denom_t = exp_cos.T.sum(dim=1) + 1e-8
        loss_ts = -torch.log((C_ts * exp_cos.T).sum(dim=1) / denom_t).mean()
        return (loss_st + loss_ts) / 2.0

    def _entropy_loss(self, R):
        log_R = torch.log(R + 1e-8)
        return self.lambda_entropy * (-(R * log_R).sum(dim=1).mean())

    def _balance_loss(self, R):
        return self.lambda_balance * R.sum(dim=0).var()

    def _cst_entropy_loss(self, C_st):
        return self.lambda_cst_entropy * (C_st * torch.log(C_st + 1e-8)).sum()

    def _mu_alignment_loss(self, mu_s, mu_t, C_st):
        return self.lambda_mu_align * (mu_s.T @ C_st - mu_t.T).norm(p='fro') ** 2


# ============================================================
# PREPROCESSING (original data.py)
# ============================================================
def preprocess_st_adata(adata, n_hvg=2000, n_pca=50):
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_cells(adata, min_genes=100)
    if n_hvg is not None and adata.n_vars > n_hvg:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor='seurat_v3')
        adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=n_pca)
    return adata

def build_spatial_weighted_graph(coords, expression, k=10):
    N = coords.shape[0]
    adj_knn = kneighbors_graph(coords, k, mode='distance', include_self=False)
    dist_matrix = adj_knn.toarray()
    nonzero = dist_matrix[dist_matrix > 0]
    sigma_s = np.median(nonzero) if len(nonzero) > 0 else 1.0
    spatial_weight = np.exp(-dist_matrix**2 / (2 * sigma_s**2))
    spatial_weight[dist_matrix == 0] = 0
    rows, cols = np.nonzero(spatial_weight)
    expr_sim = np.zeros((N, N))
    for i, j in zip(rows, cols):
        if i < j:
            try:
                corr, _ = spearmanr(expression[i], expression[j])
            except:
                corr = 0.0
            corr = np.clip(corr, -1.0, 1.0)
            if np.isnan(corr):
                corr = 0.0
            expr_sim[i, j] = expr_sim[j, i] = 0.5 * (1 + corr)
    A = spatial_weight * expr_sim
    A = np.maximum(A, 0)
    A = 0.5 * (A + A.T)
    return A

def compute_norm_adj(A):
    N = A.shape[0]
    A_tilde = A + np.eye(N)
    degrees = np.array(A_tilde.sum(axis=1)).flatten()
    d_inv_sqrt = np.where(degrees > 0, 1.0 / np.sqrt(degrees), 0.0)
    return np.diag(d_inv_sqrt) @ A_tilde @ np.diag(d_inv_sqrt)

def prepare_single_slice(adata, k=10):
    coords = adata.obsm['spatial']
    expression = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    X = adata.obsm['X_pca'].astype(np.float32)
    A = build_spatial_weighted_graph(coords, expression, k=k)
    degrees = np.sum(A, axis=1)
    total_edge_weight = 0.5 * np.sum(degrees)
    A_norm = compute_norm_adj(A)
    B = A - np.outer(degrees, degrees) / (2 * total_edge_weight + 1e-12)
    return {
        'X': torch.tensor(X, dtype=torch.float32),
        'A': A,
        'A_tensor': torch.tensor(A, dtype=torch.float32),
        'A_norm': torch.tensor(A_norm, dtype=torch.float32),
        'B': torch.tensor(B, dtype=torch.float32),
        'degrees': torch.tensor(degrees, dtype=torch.float32),
        'total_edge_weight': float(total_edge_weight),
        'N': X.shape[0],
    }

def generate_anchors_mnn(data_s, data_t, k_nn=10, top_pct=0.05):
    X_s = F.normalize(data_s['X'], p=2, dim=1)
    X_t = F.normalize(data_t['X'], p=2, dim=1)
    cos_sim = X_s @ X_t.T
    N_s, N_t = cos_sim.shape
    n_anchors = max(50, int(N_s * top_pct))
    _, top_idx_st = torch.topk(cos_sim, k=k_nn, dim=1)
    _, top_idx_ts = torch.topk(cos_sim.t(), k=k_nn, dim=1)
    mnn_pairs = []
    for i in range(N_s):
        for j_idx in range(k_nn):
            j = top_idx_st[i, j_idx].item()
            if i in top_idx_ts[j].tolist():
                mnn_pairs.append((i, j, cos_sim[i, j].item()))
    if len(mnn_pairs) == 0:
        flat_idx = np.argpartition(-cos_sim.flatten(), k_nn - 1)[:k_nn]
        for idx in flat_idx:
            i, j = idx // N_t, idx % N_t
            mnn_pairs.append((i, j, cos_sim[i, j].item()))
    mnn_pairs.sort(key=lambda x: x[2], reverse=True)
    mnn_pairs = mnn_pairs[:n_anchors * 3]
    return torch.tensor([(i, j) for i, j, _ in mnn_pairs], dtype=torch.long)


# ============================================================
# TRAINING (exact replica of train_model_a)
# ============================================================
def train_model_a(data_s, data_t, anchor_pairs, n_comm, seed, input_dim=None):
    if input_dim is None:
        input_dim = data_s['X'].shape[1]
    torch.manual_seed(seed); np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    model = CHNA_ST_Model_Original(
        input_dim=input_dim, hidden_dim=64, output_dim=32,
        n_communities_s=n_comm, n_communities_t=n_comm,
        temperature=CONFIG['temperature'],
        lambda_entropy=CONFIG['lambda_entropy'],
        lambda_balance=CONFIG['lambda_balance'],
        lambda_cst_entropy=CONFIG['lambda_cst_entropy'],
        lambda_mu_align=CONFIG['lambda_mu_align'],
    ).to(DEVICE)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=CONFIG['n_epochs'])

    A_norm_s = data_s["A_norm"].to(DEVICE)
    X_s = data_s["X"].to(DEVICE)
    B_s = data_s["B"].to(DEVICE)
    A_norm_t = data_t["A_norm"].to(DEVICE)
    X_t = data_t["X"].to(DEVICE)
    B_t = data_t["B"].to(DEVICE)
    A_s_d = data_s["A_tensor"].to(DEVICE)
    A_t_d = data_t["A_tensor"].to(DEVICE)
    anchor = anchor_pairs.to(DEVICE)

    for epoch in range(1, CONFIG['n_epochs'] + 1):
        model.train()
        optimizer.zero_grad()
        outputs = model.forward_dual_slice(
            A_norm_s, X_s, B_s, data_s["total_edge_weight"],
            A_norm_t, X_t, B_t, data_t["total_edge_weight"], anchor)

        loss_domain = CONFIG['lambda_domain'] * (outputs["loss_domain_s"] + outputs["loss_domain_t"])
        loss_anchor = CONFIG['lambda_anchor'] * outputs["loss_anchor"]
        loss_community = CONFIG['lambda_community'] * outputs["loss_community"]
        loss_entropy = outputs["loss_entropy_s"] + outputs["loss_entropy_t"]
        loss_balance = outputs["loss_balance_s"] + outputs["loss_balance_t"]
        loss_cst_entropy = outputs["loss_cst_entropy"]
        loss_mu_align = outputs["loss_mu_align"]

        conf_s = outputs["R_s"].max(dim=1).values.detach()
        conf_t = outputs["R_t"].max(dim=1).values.detach()
        mask_s = (conf_s > 0.7).float().unsqueeze(1)
        mask_t = (conf_t > 0.7).float().unsqueeze(1)
        smooth_s = (A_s_d @ (outputs["R_s"] * mask_s) - outputs["R_s"] * mask_s).norm(p=2) ** 2
        smooth_t = (A_t_d @ (outputs["R_t"] * mask_t) - outputs["R_t"] * mask_t).norm(p=2) ** 2
        loss_laplacian = CONFIG['lambda_laplacian'] * (smooth_s + smooth_t)

        total_loss = (loss_domain + loss_anchor + loss_community
                      + loss_entropy + loss_balance
                      + loss_cst_entropy + loss_mu_align + loss_laplacian)
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        scheduler.step()

    model.eval()
    with torch.no_grad():
        outputs = model.forward_dual_slice(
            A_norm_s, X_s, B_s, data_s["total_edge_weight"],
            A_norm_t, X_t, B_t, data_t["total_edge_weight"], anchor)
    return {
        'H_s': outputs['H_s'].cpu().numpy(),
        'H_t': outputs['H_t'].cpu().numpy(),
    }


# ============================================================
# EVALUATION METRICS
# ============================================================
def compute_lt(H_s, H_t, labels_s, labels_t):
    H_s_n = normalize(H_s, norm='l2', axis=1)
    H_t_n = normalize(H_t, norm='l2', axis=1)
    nn_idx = (H_t_n @ H_s_n.T).argmax(axis=1)
    transferred = labels_s[nn_idx]
    common = np.intersect1d(np.unique(labels_s), np.unique(labels_t))
    mask = np.isin(labels_t, common)
    if mask.sum() == 0:
        return 0.0
    return np.mean(labels_t[mask] == transferred[mask])

def compute_mdiag(H_s, H_t, labels_s, labels_t):
    H_s_n = normalize(H_s, norm='l2', axis=1)
    H_t_n = normalize(H_t, norm='l2', axis=1)
    common = sorted(set(labels_s) & set(labels_t))
    if len(common) < 2:
        return 0.0
    M = np.zeros((len(common), len(common)))
    for i, li in enumerate(common):
        c_s = H_s_n[labels_s == li].mean(axis=0)
        for j, lj in enumerate(common):
            c_t = H_t_n[labels_t == lj].mean(axis=0)
            M[i, j] = np.dot(c_s, c_t) / (np.linalg.norm(c_s) * np.linalg.norm(c_t) + 1e-8)
    return np.diag(M).mean()

def stratified_probe(H, labels, seed):
    idx = np.arange(len(labels))
    train_idx, test_idx = train_test_split(
        idx, test_size=0.2, random_state=seed, stratify=labels)
    clf = LogisticRegression(max_iter=2000, C=1.0, random_state=seed)
    clf.fit(H[train_idx], labels[train_idx])
    return clf.score(H[test_idx], labels[test_idx])

def zca_whiten(H):
    H_c = H - H.mean(axis=0, keepdims=True)
    cov = H_c.T @ H_c / (H.shape[0] - 1)
    U, S, _ = np.linalg.svd(cov, full_matrices=False)
    return H_c @ U @ np.diag(1.0 / np.sqrt(S + 1e-8)) @ U.T

def lga_harmonize(H_s, labels_s, H_t, labels_t):
    H_s_w = zca_whiten(H_s); H_t_w = zca_whiten(H_t)
    ls_m = np.array([LAYER_MAP.get(l, l) for l in labels_s])
    lt_m = np.array([LAYER_MAP.get(l, l) for l in labels_t])
    present = sorted(set(ls_m) & set(lt_m))
    ordered = [l for l in LAYER_ORDER if l in present]
    C_s = np.array([H_s_w[ls_m == l].mean(axis=0) for l in ordered])
    C_t = np.array([H_t_w[lt_m == l].mean(axis=0) for l in ordered])
    U, _, Vt = np.linalg.svd(C_s.T @ C_t)
    return H_s_w @ (U @ Vt), H_t_w

def random_orthogonal_matrix(d, rng):
    """Generate a random orthogonal matrix of size d×d."""
    A = rng.randn(d, d)
    Q, _ = np.linalg.qr(A)
    return Q


# ============================================================
if __name__ == "__main__":
    # MAIN EXPERIMENT
    # ============================================================
    print("\n[1] Loading data...", flush=True)
    section_ids = ["151673", "151674", "151675", "151676"]
    pairs = [("151673", "151674"), ("151674", "151675"), ("151675", "151676")]

    adatas_raw = {}
    for sid in section_ids:
        adatas_raw[sid] = sc.read_h5ad(f"{DATA_DIR}/dlpfc_{sid}.h5ad")
        print(f"  Loaded {sid}: {adatas_raw[sid].shape[0]} spots", flush=True)

    all_results = []

    for s_id, t_id in pairs:
        print(f"\n{'='*70}")
        print(f"Pair: {s_id} -> {t_id}")
        print(f"{'='*70}", flush=True)

        adata_s = preprocess_st_adata(adatas_raw[s_id].copy(),
            n_hvg=CONFIG['n_hvg'], n_pca=CONFIG['n_pca'])
        adata_t = preprocess_st_adata(adatas_raw[t_id].copy(),
            n_hvg=CONFIG['n_hvg'], n_pca=CONFIG['n_pca'])

        data_s = prepare_single_slice(adata_s, k=CONFIG['spatial_k'])
        data_t = prepare_single_slice(adata_t, k=CONFIG['spatial_k'])
        anchor_pairs = generate_anchors_mnn(data_s, data_t,
            k_nn=CONFIG['mnn_k'], top_pct=CONFIG['mnn_top_pct'])

        labels_s = np.array(adata_s.obs['layer_guess_reordered'])
        labels_t = np.array(adata_t.obs['layer_guess_reordered'])
        n_comm = CONFIG['n_communities']

        for seed in CONFIG['seeds']:
            t0 = time.time()
            print(f"  Training seed {seed}...", end=" ", flush=True)
            outputs = train_model_a(data_s, data_t, anchor_pairs, n_comm, seed)
            elapsed = time.time() - t0
            print(f"({elapsed:.0f}s)", flush=True)

            H_s = outputs['H_s']; H_t = outputs['H_t']
            D = H_s.shape[1]

            # Baseline (no rotation)
            probe_orig = stratified_probe(H_s, labels_s, seed)
            lt_orig = compute_lt(H_s, H_t, labels_s, labels_t)
            mdiag_orig = compute_mdiag(H_s, H_t, labels_s, labels_t)

            # LGA on original
            H_s_lga, H_t_w = lga_harmonize(H_s, labels_s, H_t, labels_t)
            lt_lga_orig = compute_lt(H_s_lga, H_t_w, labels_s, labels_t)
            mdiag_lga_orig = compute_mdiag(H_s_lga, H_t_w, labels_s, labels_t)

            print(f"    Baseline: Probe={probe_orig:.3f} LT={lt_orig:.3f}→{lt_lga_orig:.3f} "
                  f"Mdiag={mdiag_orig:.3f}→{mdiag_lga_orig:.3f}", flush=True)

            # Rotation experiment
            rot_results = []
            rng = np.random.RandomState(seed * 1000)  # different seed from training

            for rot_idx in range(CONFIG['n_rotations']):
                R = random_orthogonal_matrix(D, rng)
                H_s_rot = H_s @ R

                # After rotation: Probe should be PRESERVED (orthogonal transform)
                probe_rot = stratified_probe(H_s_rot, labels_s, seed)

                # After rotation: LT and Mdiag should DROP
                lt_rot = compute_lt(H_s_rot, H_t, labels_s, labels_t)
                mdiag_rot = compute_mdiag(H_s_rot, H_t, labels_s, labels_t)

                # LGA recovery
                H_s_rot_lga, _ = lga_harmonize(H_s_rot, labels_s, H_t, labels_t)
                lt_rot_lga = compute_lt(H_s_rot_lga, H_t_w, labels_s, labels_t)
                mdiag_rot_lga = compute_mdiag(H_s_rot_lga, H_t_w, labels_s, labels_t)

                rot_results.append({
                    'rot_idx': rot_idx,
                    'probe': round(float(probe_rot), 4),
                    'lt': round(float(lt_rot), 4),
                    'lt_lga': round(float(lt_rot_lga), 4),
                    'mdiag': round(float(mdiag_rot), 4),
                    'mdiag_lga': round(float(mdiag_rot_lga), 4),
                })

            # Aggregate rotation results
            lt_rot_mean = np.mean([r['lt'] for r in rot_results])
            lt_rot_std = np.std([r['lt'] for r in rot_results])
            lt_rot_lga_mean = np.mean([r['lt_lga'] for r in rot_results])
            mdiag_rot_mean = np.mean([r['mdiag'] for r in rot_results])
            mdiag_rot_std = np.std([r['mdiag'] for r in rot_results])
            mdiag_rot_lga_mean = np.mean([r['mdiag_lga'] for r in rot_results])
            probe_rot_mean = np.mean([r['probe'] for r in rot_results])

            print(f"    Rotation (avg {CONFIG['n_rotations']}×): "
                  f"Probe={probe_rot_mean:.3f} "
                  f"LT={lt_rot_mean:.3f}±{lt_rot_std:.3f}→{lt_rot_lga_mean:.3f} "
                  f"Mdiag={mdiag_rot_mean:.3f}±{mdiag_rot_std:.3f}→{mdiag_rot_lga_mean:.3f}", flush=True)

            all_results.append({
                'pair': f"{s_id}->{t_id}", 'seed': seed,
                'probe_orig': round(float(probe_orig), 4),
                'lt_orig': round(float(lt_orig), 4),
                'lt_lga_orig': round(float(lt_lga_orig), 4),
                'mdiag_orig': round(float(mdiag_orig), 4),
                'mdiag_lga_orig': round(float(mdiag_lga_orig), 4),
                'probe_rot': round(float(probe_rot_mean), 4),
                'lt_rot': round(float(lt_rot_mean), 4),
                'lt_rot_std': round(float(lt_rot_std), 4),
                'lt_rot_lga': round(float(lt_rot_lga_mean), 4),
                'mdiag_rot': round(float(mdiag_rot_mean), 4),
                'mdiag_rot_std': round(float(mdiag_rot_std), 4),
                'mdiag_rot_lga': round(float(mdiag_rot_lga_mean), 4),
                'rotations': rot_results,
                'time_s': round(float(elapsed), 1),
            })
            del outputs; gc.collect(); torch.cuda.empty_cache()


    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    print("\n" + "=" * 80)
    print("CONTROLLED ROTATION EXPERIMENT — RESULTS")
    print("=" * 80)

    print(f"\n{'Pair':<20} {'Seed':>5} {'Probe':>7} {'LT_orig':>7} {'LT_rot':>8} {'ΔLT':>7} "
          f"{'LT_LGA':>7} {'Mdiag_orig':>10} {'Mdiag_rot':>10} {'ΔMdiag':>8} {'Mdiag_LGA':>10}")
    print("-" * 105)

    for r in all_results:
        dlt = r['lt_rot'] - r['lt_orig']
        dm = r['mdiag_rot'] - r['mdiag_orig']
        print(f"{r['pair']:<20} {r['seed']:>5} {r['probe_orig']:>7.3f} {r['lt_orig']:>7.3f} "
              f"{r['lt_rot']:>8.3f} {dlt:>+7.3f} {r['lt_rot_lga']:>7.3f} "
              f"{r['mdiag_orig']:>10.3f} {r['mdiag_rot']:>10.3f} {dm:>+8.3f} {r['mdiag_rot_lga']:>10.3f}")

    # Aggregate by pair
    print(f"\n{'='*80}")
    print("AGGREGATE BY PAIR (mean ± std across seeds)")
    print(f"{'='*80}")
    print(f"{'Pair':<20} {'Probe':>8} {'LT_orig':>8} {'LT_rot':>10} {'LT_LGA':>8} "
          f"{'Mdiag_orig':>10} {'Mdiag_rot':>10} {'Mdiag_LGA':>10}")
    print("-" * 90)

    for pair_key in [f"{s}->{t}" for s, t in pairs]:
        sub = [r for r in all_results if r['pair'] == pair_key]
        means = {}
        stds = {}
        for k in ['probe_orig', 'lt_orig', 'lt_rot', 'lt_rot_lga', 'mdiag_orig', 'mdiag_rot', 'mdiag_rot_lga']:
            vals = [r[k] for r in sub]
            means[k] = np.mean(vals)
            stds[k] = np.std(vals)
        print(f"{pair_key:<20} {means['probe_orig']:>7.3f}±{stds['probe_orig']:.3f} "
              f"{means['lt_orig']:>7.3f}±{stds['lt_orig']:.3f} "
              f"{means['lt_rot']:>9.3f}±{stds['lt_rot']:.3f} "
              f"{means['lt_rot_lga']:>7.3f}±{stds['lt_rot_lga']:.3f} "
              f"{means['mdiag_orig']:>10.3f}±{stds['mdiag_orig']:.3f} "
              f"{means['mdiag_rot']:>10.3f}±{stds['mdiag_rot']:.3f} "
              f"{means['mdiag_rot_lga']:>10.3f}±{stds['mdiag_rot_lga']:.3f}")

    # Key claims
    all_probe_drop = []
    all_lt_drop = []
    all_lt_recovery = []
    all_mdiag_drop = []
    for r in all_results:
        all_probe_drop.append(abs(r['probe_rot'] - r['probe_orig']))
        all_lt_drop.append(r['lt_rot'] - r['lt_orig'])
        all_lt_recovery.append(r['lt_rot_lga'] - r['lt_rot'])
        all_mdiag_drop.append(r['mdiag_rot'] - r['mdiag_orig'])

    print(f"\n{'='*80}")
    print("KEY FINDINGS")
    print(f"{'='*80}")
    print(f"  Probe change after rotation:  {np.mean(all_probe_drop):.4f} ± {np.std(all_probe_drop):.4f}  (should be ~0)")
    print(f"  LT drop after rotation:       {np.mean(all_lt_drop):.4f} ± {np.std(all_lt_drop):.4f}  (should be <0)")
    print(f"  LT recovery after LGA:        {np.mean(all_lt_recovery):.4f} ± {np.std(all_lt_recovery):.4f}  (should be >0)")
    print(f"  Mdiag drop after rotation:    {np.mean(all_mdiag_drop):.4f} ± {np.std(all_mdiag_drop):.4f}  (should be <0)")

    # Save
    with open(f"{OUT_DIR}/rotation_results.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {OUT_DIR}/rotation_results.json", flush=True)
    print("DONE.", flush=True)