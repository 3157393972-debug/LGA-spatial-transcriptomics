#!/usr/bin/env python3
"""
Complete Natural Mismatch → LGA Recovery Experiment
====================================================
Covers AB/BC/CD for all 3 DLPFC donors (9 pairs).

Experiments per pair:
  A: Independent HVG + Shared Encoder → LGA-Oracle + LGA-Auto
  B: Shared HVG + Independent Encoder → LGA-Oracle + LGA-Auto

Measures: PD(raw), LT, Mdiag, Probe before/after each LGA variant.
"""
import warnings; warnings.filterwarnings("ignore")
import sys, os, gc, json, time
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import scanpy as sc
from scipy.spatial import ConvexHull
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize

sys.path.insert(0, os.path.dirname(__file__))
import run_rotation_historical as R

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
DATA_DIR = "/root/autodl-tmp/j"
OUT_DIR = "/root/autodl-tmp/diag_results"
os.makedirs(OUT_DIR, exist_ok=True)

SEP = "=" * 70

# ============================================================
# Auto-labeling: convex hull depth axis
# ============================================================
LAYER_ORDER = ["L1", "L2", "L3", "L4", "L5", "L6", "WM"]
LAYER_NAMES = ["Layer_1", "Layer_2", "Layer_3", "Layer_4", "Layer_5", "Layer_6", "WM"]
LAYER_MAP = {"Layer_1": "L1", "Layer_2": "L2", "Layer_3": "L3",
             "Layer_4": "L4", "Layer_5": "L5", "Layer_6": "L6", "WM": "WM"}

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


# ============================================================
# LGA-Auto: annotation-free using depth-based auto-labels
# ============================================================
def lga_auto(H_s, auto_s, H_t, auto_t):
    """LGA using auto-generated labels from spatial depth."""
    return R.lga_harmonize(H_s, auto_s, H_t, auto_t)


# ============================================================
# Independent Encoder Model
# ============================================================
class CHNA_ST_Model_Independent(nn.Module):
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
        self.encoder_s = R.SpatialEncoder(input_dim, hidden_dim, output_dim)
        self.encoder_t = R.SpatialEncoder(input_dim, hidden_dim, output_dim)
        self.soft_assign = R.SoftAssignment(temperature=temperature)
        self.center_update = R.CommunityCenter()
        self.domain_loss_fn = R.ModularityLoss()
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
        C_st = self.get_C_st(); C_ts = self.get_C_ts()
        loss_community = self._community_contrastive_loss(mu_s_new, mu_t_new, C_st, C_ts)
        return {
            "H_s": H_s, "H_t": H_t, "R_s": R_s, "R_t": R_t,
            "mu_s": mu_s_new, "mu_t": mu_t_new,
            "loss_domain_s": loss_domain_s, "loss_domain_t": loss_domain_t,
            "loss_anchor": loss_anchor, "loss_community": loss_community,
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
        cos = mu_s_n @ mu_t_n.T; exp_cos = torch.exp(cos)
        denom = exp_cos.sum(dim=1) + 1e-8
        loss_st = -torch.log((C_st * exp_cos).sum(dim=1) / denom).mean()
        denom_t = exp_cos.T.sum(dim=1) + 1e-8
        loss_ts = -torch.log((C_ts * exp_cos.T).sum(dim=1) / denom_t).mean()
        return (loss_st + loss_ts) / 2.0
    def _entropy_loss(self, R):
        return self.lambda_entropy * (-(R * torch.log(R + 1e-8)).sum(dim=1).mean())
    def _balance_loss(self, R):
        return self.lambda_balance * R.sum(dim=0).var()
    def _cst_entropy_loss(self, C_st):
        return self.lambda_cst_entropy * (C_st * torch.log(C_st + 1e-8)).sum()
    def _mu_alignment_loss(self, mu_s, mu_t, C_st):
        return self.lambda_mu_align * (mu_s.T @ C_st - mu_t.T).norm(p="fro") ** 2


# ============================================================
# Training (exact replica of train_model_a)
# ============================================================
def train_model_with_instance(model, data_s, data_t, anchor_pairs, seed):
    torch.manual_seed(seed); np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    model = model.to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=R.CONFIG["n_epochs"])
    A_norm_s = data_s["A_norm"].to(DEVICE); X_s = data_s["X"].to(DEVICE)
    B_s = data_s["B"].to(DEVICE)
    A_norm_t = data_t["A_norm"].to(DEVICE); X_t = data_t["X"].to(DEVICE)
    B_t = data_t["B"].to(DEVICE)
    A_s_d = data_s["A_tensor"].to(DEVICE); A_t_d = data_t["A_tensor"].to(DEVICE)
    anchor = anchor_pairs.to(DEVICE)
    for epoch in range(1, R.CONFIG["n_epochs"] + 1):
        model.train(); optimizer.zero_grad()
        outputs = model.forward_dual_slice(
            A_norm_s, X_s, B_s, data_s["total_edge_weight"],
            A_norm_t, X_t, B_t, data_t["total_edge_weight"], anchor)
        loss_domain = R.CONFIG["lambda_domain"] * (outputs["loss_domain_s"] + outputs["loss_domain_t"])
        loss_anchor = R.CONFIG["lambda_anchor"] * outputs["loss_anchor"]
        loss_community = R.CONFIG["lambda_community"] * outputs["loss_community"]
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
        loss_laplacian = R.CONFIG["lambda_laplacian"] * (smooth_s + smooth_t)
        total_loss = (loss_domain + loss_anchor + loss_community
                      + loss_entropy + loss_balance
                      + loss_cst_entropy + loss_mu_align + loss_laplacian)
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step(); scheduler.step()
    model.eval()
    with torch.no_grad():
        outputs = model.forward_dual_slice(
            A_norm_s, X_s, B_s, data_s["total_edge_weight"],
            A_norm_t, X_t, B_t, data_t["total_edge_weight"], anchor)
    return {"H_s": outputs["H_s"].cpu().numpy(), "H_t": outputs["H_t"].cpu().numpy()}


# ============================================================
# Metrics
# ============================================================
def procrustes_distance(H1, H2):
    H1_c = H1 - H1.mean(axis=0); H2_c = H2 - H2.mean(axis=0)
    n1, n2 = H1_c.shape[0], H2_c.shape[0]
    if n1 != n2:
        n_min = min(n1, n2)
        rng = np.random.RandomState(42)
        idx1 = rng.choice(n1, n_min, replace=False)
        idx2 = rng.choice(n2, n_min, replace=False)
        H1_c = H1_c[idx1]; H2_c = H2_c[idx2]
    U, s, Vt = np.linalg.svd(H1_c.T @ H2_c)
    return float(np.linalg.norm(H1_c @ (U @ Vt) - H2_c) / np.sqrt(H1_c.shape[0]))

def compute_lt(H_s, H_t, labels_s, labels_t):
    H_s_n = normalize(H_s, norm="l2", axis=1)
    H_t_n = normalize(H_t, norm="l2", axis=1)
    nn_idx = (H_t_n @ H_s_n.T).argmax(axis=1)
    transferred = labels_s[nn_idx]
    common = np.intersect1d(np.unique(labels_t), np.unique(transferred))
    mask = np.isin(labels_t, common) & np.isin(transferred, common)
    if mask.sum() == 0: return 0.0
    return float(np.mean(labels_t[mask] == transferred[mask]))

def compute_mdiag(H_s, H_t, labels_s, labels_t):
    H_s_n = normalize(H_s, norm="l2", axis=1)
    H_t_n = normalize(H_t, norm="l2", axis=1)
    common = sorted(set(labels_s) & set(labels_t))
    if len(common) < 2: return 0.0
    M = np.zeros((len(common), len(common)))
    for i, li in enumerate(common):
        c_s = H_s_n[labels_s == li].mean(axis=0)
        for j, lj in enumerate(common):
            c_t = H_t_n[labels_t == lj].mean(axis=0)
            M[i, j] = np.dot(c_s, c_t) / (np.linalg.norm(c_s) * np.linalg.norm(c_t) + 1e-8)
    return float(np.diag(M).mean())

def stratified_probe(H, labels, seed):
    idx = np.arange(len(labels))
    train_idx, test_idx = train_test_split(idx, test_size=0.2, random_state=seed, stratify=labels)
    clf = LogisticRegression(max_iter=2000, C=1.0, random_state=seed)
    clf.fit(H[train_idx], labels[train_idx])
    return float(clf.score(H[test_idx], labels[test_idx]))


# ============================================================
# Main experiment
# ============================================================
def run_final_experiment():
    donors = {
        "Donor1": {
            "sections": ["151507", "151508", "151509", "151510"],
            "pairs": {"AB": ("151507", "151508"), "BC": ("151508", "151509"), "CD": ("151509", "151510")},
        },
        "Donor2": {
            "sections": ["151669", "151670", "151671", "151672"],
            "pairs": {"AB": ("151669", "151670"), "BC": ("151670", "151671"), "CD": ("151671", "151672")},
        },
        "Donor3": {
            "sections": ["151673", "151674", "151675", "151676"],
            "pairs": {"AB": ("151673", "151674"), "BC": ("151674", "151675"), "CD": ("151675", "151676")},
        },
    }

    seeds = [42, 123, 456]
    all_results = []

    for donor_name, donor_info in donors.items():
        section_ids = donor_info["sections"]
        pairs = donor_info["pairs"]
        print("\n" + SEP)
        print("DONOR: {} | Sections: {}".format(donor_name, section_ids))
        print(SEP)

        adatas = {}
        for sid in section_ids:
            a = sc.read_h5ad("{}/dlpfc_{}.h5ad".format(DATA_DIR, sid))
            mask = ~a.obs["layer_guess_reordered"].isna()
            a = a[mask.values].copy() if hasattr(mask, "values") else a[mask].copy()
            adatas[sid] = a
            print("  Loaded {}: {} spots".format(sid, a.shape[0]))

        # Fit depth axis from first section
        first_sid = section_ids[0]
        first_coords = adatas[first_sid].obsm["spatial"]
        axis, intercept = fit_depth_hull_normal(first_coords)
        print("  Depth axis: ({:.4f}, {:.4f})".format(axis[0], axis[1]))

        for pair_name, (s_id, t_id) in pairs.items():
            print("\n  Pair: {} ({} -> {})".format(pair_name, s_id, t_id))
            a_s = adatas[s_id]; a_t = adatas[t_id]

            for seed in seeds:
                t0 = time.time()

                # ============================================
                # EXPERIMENT A: Independent HVG + Shared Encoder
                # ============================================
                print("    Seed {}: [A] Independent HVG...".format(seed), end=" ", flush=True)

                pp_s = R.preprocess_st_adata(a_s.copy(), n_hvg=2000, n_pca=50)
                pp_t = R.preprocess_st_adata(a_t.copy(), n_hvg=2000, n_pca=50)
                gt_s = np.array(pp_s.obs["layer_guess_reordered"])
                gt_t = np.array(pp_t.obs["layer_guess_reordered"])

                # Auto-labels for LGA-Auto
                auto_s, _ = auto_label(pp_s.obsm["spatial"], axis, intercept)
                auto_t, _ = auto_label(pp_t.obsm["spatial"], axis, intercept)

                data_s = R.prepare_single_slice(pp_s, k=10)
                data_t = R.prepare_single_slice(pp_t, k=10)
                anchors = R.generate_anchors_mnn(data_s, data_t, k_nn=10, top_pct=0.05)

                model_shared = R.CHNA_ST_Model_Original(
                    input_dim=50, hidden_dim=64, output_dim=32, n_communities_s=7,
                    temperature=R.CONFIG["temperature"], lambda_entropy=R.CONFIG["lambda_entropy"],
                    lambda_balance=R.CONFIG["lambda_balance"], lambda_cst_entropy=R.CONFIG["lambda_cst_entropy"],
                    lambda_mu_align=R.CONFIG["lambda_mu_align"])
                emb = train_model_with_instance(model_shared, data_s, data_t, anchors, seed)
                H_s_raw, H_t_raw = emb["H_s"], emb["H_t"]

                pd_raw = procrustes_distance(H_s_raw, H_t_raw)
                lt_before = compute_lt(H_s_raw, H_t_raw, gt_s, gt_t)
                md_before = compute_mdiag(H_s_raw, H_t_raw, gt_s, gt_t)
                probe_before = stratified_probe(H_s_raw, gt_s, seed)

                # LGA-Oracle
                H_s_ora, H_t_ora = R.lga_harmonize(H_s_raw, gt_s, H_t_raw, gt_t)
                lt_ora = compute_lt(H_s_ora, H_t_ora, gt_s, gt_t)
                md_ora = compute_mdiag(H_s_ora, H_t_ora, gt_s, gt_t)
                probe_ora = stratified_probe(H_s_ora, gt_s, seed)

                # LGA-Auto
                H_s_auto, H_t_auto = lga_auto(H_s_raw, auto_s, H_t_raw, auto_t)
                lt_auto = compute_lt(H_s_auto, H_t_auto, gt_s, gt_t)
                md_auto = compute_mdiag(H_s_auto, H_t_auto, gt_s, gt_t)
                probe_auto = stratified_probe(H_s_auto, gt_s, seed)

                print("PD(raw)={:.3f} LT: {:.3f}->{:.3f}(O)/{:.3f}(A) Mdiag: {:.3f}->{:.3f}(O)/{:.3f}(A) ({:.0f}s)".format(
                    pd_raw, lt_before, lt_ora, lt_auto, md_before, md_ora, md_auto, time.time() - t0))

                all_results.append({
                    "experiment": "A_independent_HVG",
                    "donor": donor_name, "pair": pair_name,
                    "source": s_id, "target": t_id, "seed": seed,
                    "pd_raw": pd_raw,
                    "before": {"LT": lt_before, "Mdiag": md_before, "Probe": probe_before},
                    "lga_oracle": {"LT": lt_ora, "Mdiag": md_ora, "Probe": probe_ora},
                    "lga_auto": {"LT": lt_auto, "Mdiag": md_auto, "Probe": probe_auto},
                })

                # ============================================
                # EXPERIMENT B: Shared HVG + Independent Encoder
                # ============================================
                t0b = time.time()
                print("    Seed {}: [B] Independent Encoder...".format(seed), end=" ", flush=True)

                combined = a_s.concatenate(a_t, batch_key="_batch")
                sc.pp.highly_variable_genes(combined, n_top_genes=2000, flavor="seurat_v3")
                shared_genes = combined.var_names[combined.var.highly_variable]
                a_s_sub = a_s[:, shared_genes].copy(); a_t_sub = a_t[:, shared_genes].copy()
                pp_s_sh = R.preprocess_st_adata(a_s_sub, n_hvg=None, n_pca=50)
                pp_t_sh = R.preprocess_st_adata(a_t_sub, n_hvg=None, n_pca=50)
                gt_s_sh = np.array(pp_s_sh.obs["layer_guess_reordered"])
                gt_t_sh = np.array(pp_t_sh.obs["layer_guess_reordered"])

                auto_s_sh, _ = auto_label(pp_s_sh.obsm["spatial"], axis, intercept)
                auto_t_sh, _ = auto_label(pp_t_sh.obsm["spatial"], axis, intercept)

                data_s_sh = R.prepare_single_slice(pp_s_sh, k=10)
                data_t_sh = R.prepare_single_slice(pp_t_sh, k=10)
                anchors_sh = R.generate_anchors_mnn(data_s_sh, data_t_sh, k_nn=10, top_pct=0.05)

                model_indep = CHNA_ST_Model_Independent(
                    input_dim=50, hidden_dim=64, output_dim=32, n_communities_s=7,
                    temperature=R.CONFIG["temperature"], lambda_entropy=R.CONFIG["lambda_entropy"],
                    lambda_balance=R.CONFIG["lambda_balance"], lambda_cst_entropy=R.CONFIG["lambda_cst_entropy"],
                    lambda_mu_align=R.CONFIG["lambda_mu_align"])
                emb_b = train_model_with_instance(model_indep, data_s_sh, data_t_sh, anchors_sh, seed)
                H_s_b, H_t_b = emb_b["H_s"], emb_b["H_t"]

                pd_raw_b = procrustes_distance(H_s_b, H_t_b)
                lt_bef_b = compute_lt(H_s_b, H_t_b, gt_s_sh, gt_t_sh)
                md_bef_b = compute_mdiag(H_s_b, H_t_b, gt_s_sh, gt_t_sh)
                probe_bef_b = stratified_probe(H_s_b, gt_s_sh, seed)

                # LGA-Oracle
                H_s_ora_b, H_t_ora_b = R.lga_harmonize(H_s_b, gt_s_sh, H_t_b, gt_t_sh)
                lt_ora_b = compute_lt(H_s_ora_b, H_t_ora_b, gt_s_sh, gt_t_sh)
                md_ora_b = compute_mdiag(H_s_ora_b, H_t_ora_b, gt_s_sh, gt_t_sh)
                probe_ora_b = stratified_probe(H_s_ora_b, gt_s_sh, seed)

                # LGA-Auto
                H_s_auto_b, H_t_auto_b = lga_auto(H_s_b, auto_s_sh, H_t_b, auto_t_sh)
                lt_auto_b = compute_lt(H_s_auto_b, H_t_auto_b, gt_s_sh, gt_t_sh)
                md_auto_b = compute_mdiag(H_s_auto_b, H_t_auto_b, gt_s_sh, gt_t_sh)
                probe_auto_b = stratified_probe(H_s_auto_b, gt_s_sh, seed)

                print("PD(raw)={:.3f} LT: {:.3f}->{:.3f}(O)/{:.3f}(A) Mdiag: {:.3f}->{:.3f}(O)/{:.3f}(A) ({:.0f}s)".format(
                    pd_raw_b, lt_bef_b, lt_ora_b, lt_auto_b, md_bef_b, md_ora_b, md_auto_b, time.time() - t0b))

                all_results.append({
                    "experiment": "B_independent_encoder",
                    "donor": donor_name, "pair": pair_name,
                    "source": s_id, "target": t_id, "seed": seed,
                    "pd_raw": pd_raw_b,
                    "before": {"LT": lt_bef_b, "Mdiag": md_bef_b, "Probe": probe_bef_b},
                    "lga_oracle": {"LT": lt_ora_b, "Mdiag": md_ora_b, "Probe": probe_ora_b},
                    "lga_auto": {"LT": lt_auto_b, "Mdiag": md_auto_b, "Probe": probe_auto_b},
                })

    # Save
    out_path = os.path.join(OUT_DIR, "natural_mismatch_lga_full.json")
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)

    # Summary
    print("\n" + SEP)
    print("SUMMARY")
    print(SEP)

    for exp_name in ["A_independent_HVG", "B_independent_encoder"]:
        exp_r = [r for r in all_results if r["experiment"] == exp_name]
        pd_raws = [r["pd_raw"] for r in exp_r]
        lt_bef = [r["before"]["LT"] for r in exp_r]
        lt_ora = [r["lga_oracle"]["LT"] for r in exp_r]
        lt_auto = [r["lga_auto"]["LT"] for r in exp_r]
        md_bef = [r["before"]["Mdiag"] for r in exp_r]
        md_ora = [r["lga_oracle"]["Mdiag"] for r in exp_r]
        md_auto = [r["lga_auto"]["Mdiag"] for r in exp_r]

        from scipy.stats import ttest_rel
        t_ora, p_ora = ttest_rel(lt_ora, lt_bef)
        t_auto, p_auto = ttest_rel(lt_auto, lt_bef)
        t_md_ora, p_md_ora = ttest_rel(md_ora, md_bef)
        t_md_auto, p_md_auto = ttest_rel(md_auto, md_bef)

        print("\n[{}] (n={})".format(exp_name, len(exp_r)))
        print("  PD(raw):    {:.4f} +/- {:.4f}".format(np.mean(pd_raws), np.std(pd_raws)))
        print("  LT before:  {:.4f}".format(np.mean(lt_bef)))
        print("  LT Oracle:  {:.4f} (t={:.2f}, p={:.6f})".format(np.mean(lt_ora), t_ora, p_ora))
        print("  LT Auto:    {:.4f} (t={:.2f}, p={:.6f})".format(np.mean(lt_auto), t_auto, p_auto))
        print("  Mdiag before: {:.4f}".format(np.mean(md_bef)))
        print("  Mdiag Oracle: {:.4f} (t={:.2f}, p={:.6f})".format(np.mean(md_ora), t_md_ora, p_md_ora))
        print("  Mdiag Auto:   {:.4f} (t={:.2f}, p={:.6f})".format(np.mean(md_auto), t_md_auto, p_md_auto))

        # By pair type
        for pair_type in ["AB", "BC", "CD"]:
            pr = [r for r in exp_r if r["pair"] == pair_type]
            if len(pr) == 0: continue
            lt_b = np.mean([r["before"]["LT"] for r in pr])
            lt_o = np.mean([r["lga_oracle"]["LT"] for r in pr])
            lt_a = np.mean([r["lga_auto"]["LT"] for r in pr])
            print("    {}: LT {:.3f} -> Oracle {:.3f} / Auto {:.3f}".format(pair_type, lt_b, lt_o, lt_a))

    print("\nResults saved to: {}".format(out_path))
    return all_results


if __name__ == "__main__":
    run_final_experiment()