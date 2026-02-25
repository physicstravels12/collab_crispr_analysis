#!/usr/bin/env python3
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.cache")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


def main():
    root = Path(__file__).resolve().parent
    cfg = yaml.safe_load((root / "config.yaml").read_text(encoding="utf-8"))
    if not bool(cfg.get("on_target_expression_filter", {}).get("enabled", True)):
        return

    out_sceptre = root / "sceptre"
    out_fig = root / "figures"
    out_sceptre.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)

    inp = out_sceptre / "on_target_results.csv"
    if not inp.exists():
        raise FileNotFoundError(f"Missing {inp}")

    de = pd.read_csv(inp)
    adata = ad.read_h5ad(cfg["data"]["h5ad_path"])

    control_col = cfg["analysis"]["control_column"]
    ctrl = adata.obs[control_col].astype(bool).values

    X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    Xc = X[ctrl]

    lib = np.asarray(Xc.sum(axis=1)).ravel().astype(float)
    scale = np.where(lib > 0, float(cfg["analysis"]["cpm_scale"]) / lib, 0.0)
    Xc_cpm = Xc.multiply(scale[:, None]).tocsr()

    control_detect_frac = np.asarray((Xc > 0).mean(axis=0)).ravel()
    control_mean_cpm = np.asarray(Xc_cpm.mean(axis=0)).ravel()

    gene_metrics = pd.DataFrame(
        {
            "gene": adata.var_names.astype(str),
            "control_detect_frac": control_detect_frac,
            "control_mean_cpm": control_mean_cpm,
        }
    )

    out = de.merge(gene_metrics, left_on="response_id", right_on="gene", how="left")
    threshold = float(cfg["on_target_expression_filter"]["threshold"])
    cpm_floor = float(cfg["on_target_expression_filter"].get("control_mean_cpm_floor", 0.0))

    out["passes_expression_filter"] = (
        (out["control_detect_frac"].fillna(0.0) >= threshold)
        & (out["control_mean_cpm"].fillna(0.0) >= cpm_floor)
    )
    out["expected_ko_direction"] = "negative"
    out["direction_matches_expectation"] = out["log_2_fold_change"].astype(float) < 0

    filtered = out[out["passes_expression_filter"]].copy()

    out.to_csv(out_sceptre / "on_target_results_with_expression_metrics.csv", index=False)
    filtered.to_csv(out_sceptre / "on_target_results_expressed_only.csv", index=False)

    stats = pd.DataFrame(
        {
            "metric": [
                "n_on_target_total",
                "n_on_target_expressed",
                "frac_expressed",
                "n_direction_match_total",
                "n_direction_match_expressed",
            ],
            "value": [
                len(out),
                len(filtered),
                float(len(filtered) / len(out)) if len(out) > 0 else np.nan,
                int(out["direction_matches_expectation"].sum()) if len(out) > 0 else 0,
                int(filtered["direction_matches_expectation"].sum()) if len(filtered) > 0 else 0,
            ],
        }
    )
    stats.to_csv(out_sceptre / "on_target_expression_filter_stats.csv", index=False)

    if len(out) > 0:
        plt.figure(figsize=(7, 5))
        x = out["control_detect_frac"].astype(float).fillna(0).values
        y = out["log_2_fold_change"].astype(float).values
        c = np.where(out["passes_expression_filter"].values, "#2563eb", "#9ca3af")
        plt.scatter(x, y, s=10, alpha=0.6, c=c, linewidths=0)
        plt.axvline(threshold, linestyle="--", color="#6b7280", linewidth=1)
        plt.axhline(0, linestyle="--", color="#6b7280", linewidth=1)
        plt.xlabel("Control detection fraction")
        plt.ylabel("On-target log2 fold-change")
        plt.title("Expression-aware on-target filtering")
        plt.tight_layout()
        plt.savefig(out_fig / "on_target_expression_aware_scatter.png", dpi=180)
        plt.close()


if __name__ == "__main__":
    main()
