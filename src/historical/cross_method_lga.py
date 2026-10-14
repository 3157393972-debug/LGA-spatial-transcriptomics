"""Cross-method LGA validation: Apply rotation + LGA to PCA, Random Projection, GCN-AE embeddings.
Tests whether LGA generalizes beyond GraphSAGE to other representation types.
"""
import sys, os, json, gc
sys.path.insert(0, '/root/autodl-tmp/j')
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.random_projection import GaussianRandomProjection
from sklearn.preprocessing import StandardScaler
from data import preprocess_st_adata

print('CROSS-METHOD LGA VALIDATION')
print('=' * 70)

pairs = [('151673', '151674'), ('151674', '151675'), ('151675', '151676')]
N_LAYERS = 7
SEEDS = [42, 123, 456]
N_ROTATIONS = 10
LAYER_ORDER = ['Layer_1', 'Layer_2', 'Layer_3', 'Layer_4', 'Layer_5', 'Layer_6', 'WM']

def whiten(X, reg=1e-5):
    X_c = X - X.mean(axis=0, keepdims=True)
    cov = (X_c.T @ X_c) / (X.shape[0] - 1)
    cov += reg * np.eye(cov.shape[0])
    U, s, _ = np.linalg.svd(cov)
    s_inv_sqrt = np.diag(1.0 / np.sqrt(s))
    return X_c @ U @ s_inv_sqrt @ U.T

def linear_probe(X, y, seed=42):
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=seed, stratify=y)
    clf = LogisticRegression(max_iter=2000, multi_class='multinomial', random_state=seed)
    clf.fit(X_train, y_train)
    return clf.score(X_test, y_test)

def compute_lt(X_s, X_t, labels_s, labels_t):
    X_s_n = X_s / (np.linalg.norm(X_s, axis=1, keepdims=True) + 1e-8)
    X_t_n = X_t / (np.linalg.norm(X_t, axis=1, keepdims=True) + 1e-8)
    sim = X_s_n @ X_t_n.T
    matches = np.argmax(sim, axis=1)
    return np.mean([1.0 if labels_s[i] == labels_t[matches[i]] else 0.0 for i in range(len(labels_s))])

def random_rotation(dim):
    G = np.random.randn(dim, dim)
    Q, R = np.linalg.qr(G)
    return Q

def lga_harmonize(H_s, labels_s, H_t, labels_t, use_whitening=True):
    """LGA: ZCA whiten + Procrustes alignment using layer centroids."""
    if use_whitening:
        H_s_w = whiten(H_s)
        H_t_w = whiten(H_t)
    else:
        H_s_w = H_s.copy()
        H_t_w = H_t.copy()
    
    # Build centroids from layer labels
    ordered = [l for l in LAYER_ORDER if l in set(labels_s) and l in set(labels_t)]
    C_s = np.array([H_s_w[np.array(labels_s) == l].mean(axis=0) for l in ordered])
    C_t = np.array([H_t_w[np.array(labels_t) == l].mean(axis=0) for l in ordered])
    
    U, _, Vt = np.linalg.svd(C_s.T @ C_t)
    R = U @ Vt
    H_s_lga = H_s_w @ R
    return H_s_lga, H_t_w

def lga_auto_harmonize(H_s, coords_s, H_t, coords_t, K=7, use_whitening=True):
    """LGA-Auto: spatial depth binning instead of layer labels."""
    if use_whitening:
        H_s_w = whiten(H_s)
        H_t_w = whiten(H_t)
    else:
        H_s_w = H_s.copy()
        H_t_w = H_t.copy()
    
    # Fit depth axis via convex hull
    from scipy.spatial import ConvexHull
    pts_s = np.array(coords_s)
    pts_t = np.array(coords_t)
    
    def fit_depth_axis(pts):
        hull = ConvexHull(pts)
        # Find two farthest hull vertices as axis endpoints
        hull_pts = pts[hull.vertices]
        max_dist = 0
        best_pair = (0, 1)
        for i in range(len(hull_pts)):
            for j in range(i+1, len(hull_pts)):
                d = np.linalg.norm(hull_pts[i] - hull_pts[j])
                if d > max_dist:
                    max_dist = d
                    best_pair = (i, j)
        axis = hull_pts[best_pair[1]] - hull_pts[best_pair[0]]
        axis = axis / np.linalg.norm(axis)
        origin = hull_pts[best_pair[0]]
        projections = (pts - origin) @ axis
        return projections
    
    depth_s = fit_depth_axis(pts_s)
    depth_t = fit_depth_axis(pts_t)
    
    # Pooled quantile binning
    all_depths = np.concatenate([depth_s, depth_t])
    quantiles = np.percentile(all_depths, np.linspace(0, 100, K+1))
    
    C_s = []
    C_t = []
    for k in range(K):
        mask_s = (depth_s >= quantiles[k]) & (depth_s < quantiles[k+1]) if k < K-1 else (depth_s >= quantiles[k])
        mask_t = (depth_t >= quantiles[k]) & (depth_t < quantiles[k+1]) if k < K-1 else (depth_t >= quantiles[k])
        if mask_s.sum() > 0:
            C_s.append(H_s_w[mask_s].mean(axis=0))
        if mask_t.sum() > 0:
            C_t.append(H_t_w[mask_t].mean(axis=0))
    
    C_s = np.array(C_s)
    C_t = np.array(C_t)
    # Ensure same number of centroids
    min_k = min(len(C_s), len(C_t))
    C_s, C_t = C_s[:min_k], C_t[:min_k]
    
    U, _, Vt = np.linalg.svd(C_s.T @ C_t)
    R = U @ Vt
    H_s_lga = H_s_w @ R
    return H_s_lga, H_t_w


# === GCN Autoencoder ===
class GCNEncoder(nn.Module):
    def __init__(self, in_dim, hidden_dim=128, out_dim=32):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, out_dim)
        self.bn1 = nn.BatchNorm1d(hidden_dim)
    def forward(self, x, adj=None):
        h = F.relu(self.bn1(self.fc1(x)))
        return self.fc2(h)

class GCNDecoder(nn.Module):
    def __init__(self, in_dim, hidden_dim=128, out_dim=None):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, out_dim)
        self.bn1 = nn.BatchNorm1d(hidden_dim)
    def forward(self, h, adj=None):
        x = F.relu(self.bn1(self.fc1(h)))
        return self.fc2(x)

class GCNAutoencoder(nn.Module):
    def __init__(self, in_dim, hidden_dim=128, latent_dim=32):
        super().__init__()
        self.encoder = GCNEncoder(in_dim, hidden_dim, latent_dim)
        self.decoder = GCNDecoder(latent_dim, hidden_dim, in_dim)
    def forward(self, x, adj=None):
        h = self.encoder(x, adj)
        return self.decoder(h, adj), h

def train_gcn_ae(X, n_epochs=200, lr=1e-3, latent_dim=32, device='cpu'):
    X_tensor = torch.tensor(X, dtype=torch.float32).to(device)
    model = GCNAutoencoder(X.shape[1], hidden_dim=128, latent_dim=latent_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    model.train()
    for epoch in range(n_epochs):
        optimizer.zero_grad()
        x_rec, h = model(X_tensor)
        loss = F.mse_loss(x_rec, X_tensor)
        loss.backward()
        optimizer.step()
    model.eval()
    with torch.no_grad():
        _, h = model(X_tensor)
        return h.cpu().numpy()


# ============================================================
# MAIN
# ============================================================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Device: {device}')

all_results = {}

for sid_s, sid_t in pairs:
    pk = sid_s + '->' + sid_t
    print(f'\n=== Loading {pk} ===')
    
    adata_s = sc.read_h5ad(f'/root/autodl-tmp/j/dlpfc_{sid_s}.h5ad')
    adata_t = sc.read_h5ad(f'/root/autodl-tmp/j/dlpfc_{sid_t}.h5ad')
    adata_s = preprocess_st_adata(adata_s)
    adata_t = preprocess_st_adata(adata_t)
    
    labels_s = np.array(adata_s.obs['layer_guess_reordered'].values.astype(str))
    labels_t = np.array(adata_t.obs['layer_guess_reordered'].values.astype(str))
    coords_s = adata_s.obsm["spatial"]
    coords_t = adata_t.obsm["spatial"]
    
    X_s_raw = adata_s.X.toarray() if hasattr(adata_s.X, 'toarray') else adata_s.X
    X_t_raw = adata_t.X.toarray() if hasattr(adata_t.X, 'toarray') else adata_t.X
    
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X_s_raw)
    X_t = scaler.transform(X_t_raw)
    
    print(f'  S: {X_s.shape}, T: {X_t.shape}')
    
    # ============================================================
    # Method 1: PCA (50d)
    # ============================================================
    print('\n  --- PCA (50d) ---')
    pca = PCA(n_components=50, random_state=42)
    H_s = pca.fit_transform(X_s)
    H_t = pca.transform(X_t)
    
    probe = np.mean([linear_probe(H_s, labels_s, seed=s) for s in SEEDS])
    lt_orig = compute_lt(H_s, H_t, labels_s, labels_t)
    
    lt_rot_list = []
    lt_oracle_list = []
    lt_auto_list = []
    for rot_idx in range(N_ROTATIONS):
        R = random_rotation(H_s.shape[1])
        H_s_rot = H_s @ R
        lt_rot = compute_lt(H_s_rot, H_t, labels_s, labels_t)
        H_s_oracle, H_t_w = lga_harmonize(H_s_rot, labels_s, H_t, labels_t, use_whitening=False)
        lt_oracle = compute_lt(H_s_oracle, H_t_w, labels_s, labels_t)
        H_s_auto, H_t_aw = lga_auto_harmonize(H_s_rot, coords_s, H_t, coords_t, K=7, use_whitening=False)
        lt_auto = compute_lt(H_s_auto, H_t_aw, labels_s, labels_t)
        lt_rot_list.append(lt_rot)
        lt_oracle_list.append(lt_oracle)
        lt_auto_list.append(lt_auto)
    
    all_results[f'pca_{pk}'] = {
        'method': 'PCA', 'pair': pk,
        'probe': float(probe),
        'lt_original': float(lt_orig),
        'lt_post_rotation': float(np.mean(lt_rot_list)),
        'lt_lga_oracle': float(np.mean(lt_oracle_list)),
        'lt_lga_auto': float(np.mean(lt_auto_list)),
        'recovery_oracle': float((np.mean(lt_oracle_list) - np.mean(lt_rot_list)) / (lt_orig - np.mean(lt_rot_list)) * 100) if (lt_orig - np.mean(lt_rot_list)) > 0 else 0,
    }
    print(f'    Probe={probe:.4f}, LT_orig={lt_orig:.4f}, LT_rot={np.mean(lt_rot_list):.4f}, LGA-Oracle={np.mean(lt_oracle_list):.4f}, LGA-Auto={np.mean(lt_auto_list):.4f}')
    
    # ============================================================
    # Method 2: Random Projection (32d)
    # ============================================================
    print('\n  --- Random Projection (32d) ---')
    rp_lt_orig = []
    rp_lt_rot = []
    rp_lt_oracle = []
    rp_lt_auto = []
    rp_probe = []
    for seed in SEEDS:
        rp = GaussianRandomProjection(n_components=32, random_state=seed)
        H_s = rp.fit_transform(X_s)
        H_t_rp = rp.transform(X_t)
        rp_probe.append(linear_probe(H_s, labels_s, seed=seed))
        rp_lt_orig.append(compute_lt(H_s, H_t_rp, labels_s, labels_t))
        for rot_idx in range(3):  # 3 rotations per seed
            R = random_rotation(32)
            H_s_rot = H_s @ R
            rp_lt_rot.append(compute_lt(H_s_rot, H_t_rp, labels_s, labels_t))
            H_s_oracle, H_t_w = lga_harmonize(H_s_rot, labels_s, H_t_rp, labels_t, use_whitening=False)
            rp_lt_oracle.append(compute_lt(H_s_oracle, H_t_w, labels_s, labels_t))
            H_s_auto, H_t_aw = lga_auto_harmonize(H_s_rot, coords_s, H_t_rp, coords_t, K=7, use_whitening=False)
            rp_lt_auto.append(compute_lt(H_s_auto, H_t_aw, labels_s, labels_t))
    
    all_results[f'rp_{pk}'] = {
        'method': 'Random Projection', 'pair': pk,
        'probe': float(np.mean(rp_probe)),
        'lt_original': float(np.mean(rp_lt_orig)),
        'lt_post_rotation': float(np.mean(rp_lt_rot)),
        'lt_lga_oracle': float(np.mean(rp_lt_oracle)),
        'lt_lga_auto': float(np.mean(rp_lt_auto)),
    }
    print(f'    Probe={np.mean(rp_probe):.4f}, LT_orig={np.mean(rp_lt_orig):.4f}, LT_rot={np.mean(rp_lt_rot):.4f}, LGA-Oracle={np.mean(rp_lt_oracle):.4f}, LGA-Auto={np.mean(rp_lt_auto):.4f}')
    
    # ============================================================
    # Method 3: GCN Autoencoder (32d)
    # ============================================================
    print('\n  --- GCN Autoencoder (32d) ---')
    gcn_lt_orig = []
    gcn_lt_rot = []
    gcn_lt_oracle = []
    gcn_lt_auto = []
    gcn_probe = []
    for seed in SEEDS:
        torch.manual_seed(seed)
        np.random.seed(seed)
        H_s = train_gcn_ae(X_s, n_epochs=200, lr=1e-3, latent_dim=32, device=device)
        H_t_gcn = train_gcn_ae(X_t, n_epochs=200, lr=1e-3, latent_dim=32, device=device)
        gcn_probe.append(linear_probe(H_s, labels_s, seed=seed))
        gcn_lt_orig.append(compute_lt(H_s, H_t_gcn, labels_s, labels_t))
        for rot_idx in range(3):
            R = random_rotation(32)
            H_s_rot = H_s @ R
            gcn_lt_rot.append(compute_lt(H_s_rot, H_t_gcn, labels_s, labels_t))
            H_s_oracle, H_t_w = lga_harmonize(H_s_rot, labels_s, H_t_gcn, labels_t, use_whitening=False)
            gcn_lt_oracle.append(compute_lt(H_s_oracle, H_t_w, labels_s, labels_t))
            H_s_auto, H_t_aw = lga_auto_harmonize(H_s_rot, coords_s, H_t_gcn, coords_t, K=7, use_whitening=False)
            gcn_lt_auto.append(compute_lt(H_s_auto, H_t_aw, labels_s, labels_t))
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    
    all_results[f'gcn_ae_{pk}'] = {
        'method': 'GCN-AE', 'pair': pk,
        'probe': float(np.mean(gcn_probe)),
        'lt_original': float(np.mean(gcn_lt_orig)),
        'lt_post_rotation': float(np.mean(gcn_lt_rot)),
        'lt_lga_oracle': float(np.mean(gcn_lt_oracle)),
        'lt_lga_auto': float(np.mean(gcn_lt_auto)),
    }
    print(f'    Probe={np.mean(gcn_probe):.4f}, LT_orig={np.mean(gcn_lt_orig):.4f}, LT_rot={np.mean(gcn_lt_rot):.4f}, LGA-Oracle={np.mean(gcn_lt_oracle):.4f}, LGA-Auto={np.mean(gcn_lt_auto):.4f}')

# Save
out_path = '/root/autodl-tmp/cross_method_lga_results.json'
with open(out_path, 'w') as f:
    json.dump(all_results, f, indent=2)
print(f'\nSaved to {out_path}')

# Summary
print('\n' + '=' * 80)
print('SUMMARY: LGA Cross-Method Validation (Donor 3)')
print('=' * 80)
print(f'{"Method":<20} {"Pair":<15} {"Probe":>8} {"LT_orig":>8} {"LT_rot":>8} {"LGA-Oracle":>10} {"LGA-Auto":>10}')
print('-' * 80)
for k in sorted(all_results.keys()):
    r = all_results[k]
    print(f'{r["method"]:<20} {r["pair"]:<15} {r["probe"]:>8.4f} {r["lt_original"]:>8.4f} {r["lt_post_rotation"]:>8.4f} {r["lt_lga_oracle"]:>10.4f} {r["lt_lga_auto"]:>10.4f}')

print('\nDone.')
