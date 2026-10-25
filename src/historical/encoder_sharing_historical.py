#!/usr/bin/env python3
"""
Experiment: Shared vs Independent Encoder
=========================================
Tests whether the encoder is truly shared, and what happens when we make it independent.

Design:
  1. Train with SHARED encoder (current code: encoder_t = encoder_s)
  2. Train with INDEPENDENT encoder (separate SpatialEncoder instances)
  3. Both use identical preprocessing, anchors, seeds
  4. Measure: PD, CD, LT

Hypothesis: Independent encoder produces larger PD,
  demonstrating that the current "shared" design actively suppresses
  coordinate incompatibility.
"""
import warnings; warnings.filterwarnings("ignore")
import sys, os, gc, json, time
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity

sys.path.insert(0, os.path.dirname(__file__))
import run_rotation_historical as R

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
DATA_DIR = "/root/autodl-tmp/j"
OUT_DIR = "/root/autodl-tmp/diag_results"
os.makedirs(OUT_DIR, exist_ok=True)

SEP = "=" * 70

# ============================================================
# Independent Encoder Model (encoder_t != encoder_s)
# ============================================================
class CHNA_ST_Model_Independent(nn.Module):
    """Same architecture, but SEPARATE encoder instances for source and target."""
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

        # KEY CHANGE: separate encoder instances
        self.encoder_s = R.SpatialEncoder(input_dim, hidden_dim, output_dim)
        self.encoder_t = R.SpatialEncoder(input_dim, hidden_dim, output_dim)  # INDEPENDENT!
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
        return self.lambda_mu_align * (mu_s.T @ C_st - mu_t.T).norm(p="fro") ** 2


# ============================================================
# Training function — EXACT replica of train_model_a
# but accepts a model instance so we can swap shared vs independent
# ============================================================
def train_model_with_instance(model, data_s, data_t, anchor_pairs, seed):
    """Exact replica of run_rotation.train_model_a, but with model as parameter."""
    torch.manual_seed(seed); np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    model = model.to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=R.CONFIG["n_epochs"])

    A_norm_s = data_s["A_norm"].to(DEVICE)
    X_s = data_s["X"].to(DEVICE)
    B_s = data_s["B"].to(DEVICE)
    A_norm_t = data_t["A_norm"].to(DEVICE)
    X_t = data_t["X"].to(DEVICE)
    B_t = data_t["B"].to(DEVICE)
    A_s_d = data_s["A_tensor"].to(DEVICE)
    A_t_d = data_t["A_tensor"].to(DEVICE)
    anchor = anchor_pairs.to(DEVICE)

    for epoch in range(1, R.CONFIG["n_epochs"] + 1):
        model.train()
        optimizer.zero_grad()
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
        optimizer.step()
        scheduler.step()

    model.eval()
    with torch.no_grad():
        outputs = model.forward_dual_slice(
            A_norm_s, X_s, B_s, data_s["total_edge_weight"],
            A_norm_t, X_t, B_t, data_t["total_edge_weight"], anchor)
    return {
        "H_s": outputs["H_s"].cpu().numpy(),
        "H_t": outputs["H_t"].cpu().numpy(),
    }


# ============================================================
# Metrics
# ============================================================
def procrustes_distance(H1, H2):
    H1_c = H1 - H1.mean(axis=0)
    H2_c = H2 - H2.mean(axis=0)
    n1, n2 = H1_c.shape[0], H2_c.shape[0]
    if n1 != n2:
        n_min = min(n1, n2)
        rng = np.random.RandomState(42)
        idx1 = rng.choice(n1, n_min, replace=False)
        idx2 = rng.choice(n2, n_min, replace=False)
        H1_c = H1_c[idx1]
        H2_c = H2_c[idx2]
    U, s, Vt = np.linalg.svd(H1_c.T @ H2_c)
    R_opt = U @ Vt
    return float(np.linalg.norm(H1_c @ R_opt - H2_c) / np.sqrt(H1_c.shape[0]))

def covariance_distance(H1, H2):
    C1 = np.cov(H1.T)
    C2 = np.cov(H2.T)
    return float(np.linalg.norm(C1 - C2, "fro"))

def compute_lt(H_s, H_t, gt_s, gt_t):
    sim = cosine_similarity(H_t, H_s)
    nn = np.argmax(sim, axis=1)
    return float(np.mean(gt_t == gt_s[nn]))


# ============================================================
# Main experiment
# ============================================================
def run_encoder_sharing_experiment():
    donors = {
        "Donor1": (["151507", "151508", "151509"],
                   [("151507", "151508"), ("151508", "151509"), ("151507", "151509")]),
        "Donor2": (["151669", "151670", "151671"],
                   [("151669", "151670"), ("151670", "151671"), ("151669", "151671")]),
        "Donor3": (["151673", "151674", "151675"],
                   [("151673", "151674"), ("151674", "151675"), ("151673", "151675")]),
    }

    seeds = [42, 123, 456, 789, 1024]
    all_results = {
        "config": {
            "n_hvg": 2000, "n_pca": 50,
            "input_dim": 50, "hidden_dim": 64, "output_dim": 32,
            "n_communities": 7,
            "n_epochs": 500, "lr": 1e-3, "weight_decay": 1e-4,
            "seeds": seeds,
        },
        "results": []
    }

    for donor_name, (section_ids, pairs) in donors.items():
        print(f"\n{SEP}")
        print(f"DONOR: {donor_name} | Sections: {section_ids}")
        print(SEP)

        adatas = {}
        for sid in section_ids:
            a = sc.read_h5ad(f"{DATA_DIR}/dlpfc_{sid}.h5ad")
            mask = ~a.obs["layer_guess_reordered"].isna()
            a = a[mask.values].copy() if hasattr(mask, "values") else a[mask].copy()
            adatas[sid] = a
            print(f"  Loaded {sid}: {a.shape[0]} spots")

        for s_id, t_id in pairs:
            print(f"\n  Pair: {s_id} -> {t_id}")

            # Preprocess: SHARED HVG (concatenate -> joint HVG -> subset)
            combined = adatas[s_id].concatenate(adatas[t_id], batch_key="_batch")
            sc.pp.highly_variable_genes(combined, n_top_genes=2000, flavor="seurat_v3")
            shared_genes = combined.var_names[combined.var.highly_variable]
            a_s_sub = adatas[s_id][:, shared_genes].copy()
            a_t_sub = adatas[t_id][:, shared_genes].copy()
            adata_s = R.preprocess_st_adata(a_s_sub, n_hvg=None, n_pca=50)
            adata_t = R.preprocess_st_adata(a_t_sub, n_hvg=None, n_pca=50)

            gt_s = np.array(adata_s.obs["layer_guess_reordered"])
            gt_t = np.array(adata_t.obs["layer_guess_reordered"])

            data_s = R.prepare_single_slice(adata_s, k=10)
            data_t = R.prepare_single_slice(adata_t, k=10)
            anchor_pairs = R.generate_anchors_mnn(data_s, data_t, k_nn=10, top_pct=0.05)

            for seed in seeds:
                t0 = time.time()

                # --- SHARED ENCODER (original) ---
                print(f"    Seed {seed}: shared...", end=" ", flush=True)
                model_shared = R.CHNA_ST_Model_Original(
                    input_dim=50, hidden_dim=64, output_dim=32,
                    n_communities_s=7, temperature=R.CONFIG["temperature"],
                    lambda_entropy=R.CONFIG["lambda_entropy"],
                    lambda_balance=R.CONFIG["lambda_balance"],
                    lambda_cst_entropy=R.CONFIG["lambda_cst_entropy"],
                    lambda_mu_align=R.CONFIG["lambda_mu_align"],
                )
                emb_shared = train_model_with_instance(model_shared, data_s, data_t, anchor_pairs, seed)
                pd_shared = procrustes_distance(emb_shared["H_s"], emb_shared["H_t"])
                cd_shared = covariance_distance(emb_shared["H_s"], emb_shared["H_t"])
                lt_shared = compute_lt(emb_shared["H_s"], emb_shared["H_t"], gt_s, gt_t)
                t_shared = time.time() - t0

                # --- INDEPENDENT ENCODER ---
                print(f"independent...", end=" ", flush=True)
                model_indep = CHNA_ST_Model_Independent(
                    input_dim=50, hidden_dim=64, output_dim=32,
                    n_communities_s=7, temperature=R.CONFIG["temperature"],
                    lambda_entropy=R.CONFIG["lambda_entropy"],
                    lambda_balance=R.CONFIG["lambda_balance"],
                    lambda_cst_entropy=R.CONFIG["lambda_cst_entropy"],
                    lambda_mu_align=R.CONFIG["lambda_mu_align"],
                )
                emb_indep = train_model_with_instance(model_indep, data_s, data_t, anchor_pairs, seed)
                pd_indep = procrustes_distance(emb_indep["H_s"], emb_indep["H_t"])
                cd_indep = covariance_distance(emb_indep["H_s"], emb_indep["H_t"])
                lt_indep = compute_lt(emb_indep["H_s"], emb_indep["H_t"], gt_s, gt_t)
                t_indep = time.time() - t0 - t_shared

                result = {
                    "donor": donor_name,
                    "pair": f"{s_id}-{t_id}",
                    "seed": seed,
                    "shared": {"PD": pd_shared, "CD": cd_shared, "LT": lt_shared},
                    "independent": {"PD": pd_indep, "CD": cd_indep, "LT": lt_indep},
                }
                all_results["results"].append(result)

                print(f"(shared: {t_shared:.0f}s, indep: {t_indep:.0f}s)")
                print(f"      Shared:     PD={pd_shared:.4f}  CD={cd_shared:.4f}  LT={lt_shared:.4f}")
                print(f"      Independent: PD={pd_indep:.4f}  CD={cd_indep:.4f}  LT={lt_indep:.4f}")
                print(f"      Delta PD = {pd_indep - pd_shared:+.4f}")

    # Save results
    out_path = os.path.join(OUT_DIR, "encoder_sharing.json")
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)

    # Summary
    print(f"\n{SEP}")
    print("SUMMARY")
    print(SEP)
    results = all_results["results"]
    shared_pds = [r["shared"]["PD"] for r in results]
    indep_pds = [r["independent"]["PD"] for r in results]
    shared_lts = [r["shared"]["LT"] for r in results]
    indep_lts = [r["independent"]["LT"] for r in results]
    shared_cds = [r["shared"]["CD"] for r in results]
    indep_cds = [r["independent"]["CD"] for r in results]

    print(f"\nShared encoder (n={len(shared_pds)}):")
    print(f"  PD: {np.mean(shared_pds):.4f} +/- {np.std(shared_pds):.4f}")
    print(f"  CD: {np.mean(shared_cds):.4f} +/- {np.std(shared_cds):.4f}")
    print(f"  LT: {np.mean(shared_lts):.4f} +/- {np.std(shared_lts):.4f}")
    print(f"\nIndependent encoder (n={len(indep_pds)}):")
    print(f"  PD: {np.mean(indep_pds):.4f} +/- {np.std(indep_pds):.4f}")
    print(f"  CD: {np.mean(indep_cds):.4f} +/- {np.std(indep_cds):.4f}")
    print(f"  LT: {np.mean(indep_lts):.4f} +/- {np.std(indep_lts):.4f}")

    from scipy.stats import ttest_rel
    t_pd, p_pd = ttest_rel(indep_pds, shared_pds)
    t_lt, p_lt = ttest_rel(indep_lts, shared_lts)
    t_cd, p_cd = ttest_rel(indep_cds, shared_cds)
    print(f"\nPaired t-test (independent vs shared):")
    print(f"  PD: t={t_pd:.4f}, p={p_pd:.6f}")
    print(f"  CD: t={t_cd:.4f}, p={p_cd:.6f}")
    print(f"  LT: t={t_lt:.4f}, p={p_lt:.6f}")

    print(f"\nResults saved to: {out_path}")
    return all_results


if __name__ == "__main__":
    run_encoder_sharing_experiment()