#!/usr/bin/env python3
from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.cache")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import yaml
from statsmodels.stats.multitest import multipletests


def log1p_cpm_sparse(X_csr: sp.csr_matrix, cpm_scale: float) -> sp.csr_matrix:
    libsize = np.asarray(X_csr.sum(axis=1)).ravel().astype(np.float64)
    scale = np.where(libsize > 0, cpm_scale / libsize, 0.0)
    X_cpm = X_csr.multiply(scale[:, None]).tocsr().astype(np.float32)
    X_log = X_cpm.copy()
    X_log.data = np.log1p(X_log.data)
    return X_log


def main():
    root = Path(__file__).resolve().parent
    cfg = yaml.safe_load((root / "config.yaml").read_text(encoding="utf-8"))

    out_dir = root / "perturbation_burden"
    fig_dir = root / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(cfg["data"]["h5ad_path"])
    obs = adata.obs.copy()

    if sp.issparse(adata.X):
        X = adata.X.tocsr()
    else:
        X = sp.csr_matrix(adata.X)

    disc_path = root / "sceptre" / "trans_panel_results.csv"
    if not disc_path.exists():
        raise FileNotFoundError(f"Missing {disc_path}")
    disc = pd.read_csv(disc_path)
    if len(disc) == 0:
        alt = root / "sceptre" / "on_target_results.csv"
        if alt.exists():
            disc = pd.read_csv(alt)

    if "p_value_adjusted" not in disc.columns and "p_value" in disc.columns and len(disc) > 0:
        _, padj, _, _ = multipletests(disc["p_value"].astype(float).values, method="fdr_bh")
        disc["p_value_adjusted"] = padj

    alpha = float(cfg["analysis"]["fdr_alpha"])
    responsive = disc.loc[disc["p_value_adjusted"] < alpha, "response_id"].astype(str).unique().tolist()

    # fallback in case nothing passes: pick top 200 by p-value
    if len(responsive) == 0 and len(disc) > 0:
        responsive = (
            disc.sort_values("p_value")
            .head(200)["response_id"]
            .astype(str)
            .unique()
            .tolist()
        )

    var_names = adata.var_names.astype(str).tolist()
    gene_to_idx = {g: i for i, g in enumerate(var_names)}
    idx = [gene_to_idx[g] for g in responsive if g in gene_to_idx]

    if len(idx) == 0:
        pd.DataFrame().to_csv(out_dir / "per_cell_scores.csv", index=False)
        (out_dir / "burden_summary.txt").write_text("No responsive genes available for burden calculation.\n", encoding="utf-8")
        return

    X_log = log1p_cpm_sparse(X, float(cfg["analysis"]["cpm_scale"]))[:, idx]
    ctrl_mask = obs[cfg["analysis"]["control_column"]].astype(bool).values

    X_ctrl = X_log[ctrl_mask]
    mu = np.asarray(X_ctrl.mean(axis=0)).ravel()
    mu2 = np.asarray(X_ctrl.power(2).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(mu2 - mu ** 2, 1e-8))

    # z-score per cell x responsive gene
    n_cells = X_log.shape[0]
    Z = X_log.toarray()
    Z = (Z - mu[None, :]) / sd[None, :]

    z_abs_thresh = float(cfg["perturbation_burden"]["z_abs_threshold"])
    burden_n = (np.abs(Z) >= z_abs_thresh).sum(axis=1)
    burden_weighted = np.mean(np.abs(Z), axis=1)

    out = obs.copy()
    out["cell_id"] = out.index.astype(str)
    out["perturbation_burden_n"] = burden_n
    out["perturbation_burden_weighted"] = burden_weighted
    out["n_responsive_genes_used"] = len(idx)

    cols = [
        "cell_id",
        cfg["analysis"]["control_column"],
        "assigned_guide",
        "perturbation_burden_n",
        "perturbation_burden_weighted",
        "n_responsive_genes_used",
    ]
    for c in ["lane", "sample_id", "num_umis", "num_features"]:
        if c in out.columns:
            cols.append(c)
    out[cols].to_csv(out_dir / "per_cell_scores.csv", index=False)

    summary = (
        out.groupby(cfg["analysis"]["control_column"])[["perturbation_burden_n", "perturbation_burden_weighted"]]
        .agg(["mean", "median", "std"])
    )
    summary.to_csv(out_dir / "burden_group_summary.csv")

    plt.figure(figsize=(6, 4))
    sns.boxplot(
        data=out,
        x=cfg["analysis"]["control_column"],
        y="perturbation_burden_n",
        color="#93c5fd",
    )
    plt.title("Per-cell perturbation burden (count of |z|>=2 genes)")
    plt.tight_layout()
    plt.savefig(fig_dir / "perturbation_burden_n_boxplot.png", dpi=180)
    plt.close()

    plt.figure(figsize=(6, 4))
    sns.boxplot(
        data=out,
        x=cfg["analysis"]["control_column"],
        y="perturbation_burden_weighted",
        color="#86efac",
    )
    plt.title("Per-cell perturbation burden (mean |z|)")
    plt.tight_layout()
    plt.savefig(fig_dir / "perturbation_burden_weighted_boxplot.png", dpi=180)
    plt.close()

    (out_dir / "burden_summary.txt").write_text(
        "\n".join(
            [
                f"responsive_gene_count={len(idx)}",
                f"control_mean_burden_n={out.loc[out[cfg['analysis']['control_column']] == True, 'perturbation_burden_n'].mean():.3f}",
                f"perturbed_mean_burden_n={out.loc[out[cfg['analysis']['control_column']] == False, 'perturbation_burden_n'].mean():.3f}",
                f"control_mean_burden_weighted={out.loc[out[cfg['analysis']['control_column']] == True, 'perturbation_burden_weighted'].mean():.3f}",
                f"perturbed_mean_burden_weighted={out.loc[out[cfg['analysis']['control_column']] == False, 'perturbation_burden_weighted'].mean():.3f}",
            ]
        ) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
