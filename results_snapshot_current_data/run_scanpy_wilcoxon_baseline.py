#!/usr/bin/env python3
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.cache")
os.environ.setdefault("NUMBA_CACHE_DIR", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.cache/numba")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


def parse_target_from_guide(guide: str) -> str:
    g = str(guide)
    if "NO_SITE" in g or "NON-GENE_SITE" in g:
        return "non-targeting"
    parts = g.split("_")
    if g.startswith("Cas9_5000_") and len(parts) >= 3:
        return parts[2]
    if g.startswith("A4_") and len(parts) >= 2:
        return parts[1]
    if len(parts) >= 2:
        return parts[1]
    return g


def safe_volcano(df: pd.DataFrame, out_png: Path, title: str, alpha: float):
    if len(df) == 0:
        return
    p = np.clip(df.get("pvals_adj", pd.Series(np.ones(len(df)))).astype(float).values, 1e-300, 1.0)
    lfc = df.get("logfoldchanges", pd.Series(np.zeros(len(df)))).astype(float).values
    sig = p < alpha

    plt.figure(figsize=(7, 5))
    plt.scatter(lfc, -np.log10(p), s=5, alpha=0.35, c="#6b7280", linewidths=0)
    if np.any(sig):
        plt.scatter(lfc[sig], -np.log10(p[sig]), s=7, alpha=0.65, c="#dc2626", linewidths=0)
    plt.axhline(-np.log10(alpha), color="#9ca3af", linestyle="--", linewidth=1)
    plt.xlabel("log fold-change")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def main():
    root = Path(__file__).resolve().parent
    cfg = yaml.safe_load((root / "config.yaml").read_text(encoding="utf-8"))
    if not bool(cfg.get("scanpy_de", {}).get("enabled", True)):
        return

    out_scanpy = root / "scanpy"
    out_fig = root / "figures"
    out_scanpy.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(cfg["data"]["h5ad_path"])
    adata = adata.copy()

    control_col = cfg["analysis"]["control_column"]
    guide_col = cfg["analysis"]["assigned_guide_column"]
    alpha = float(cfg["scanpy_de"]["fdr_alpha"])
    method = str(cfg["scanpy_de"].get("method", "wilcoxon"))
    min_cells_target = int(cfg["scanpy_de"].get("min_cells_per_target", cfg["analysis"]["min_cells_per_target"]))
    top_n_export = int(cfg["scanpy_de"].get("top_n_genes_export", 300))

    adata.obs[guide_col] = adata.obs[guide_col].astype(str)
    adata.obs["parsed_target"] = adata.obs[guide_col].map(parse_target_from_guide)
    adata.obs["de_global_group"] = np.where(adata.obs[control_col].astype(bool), "control", "perturbed")

    sc.pp.normalize_total(adata, target_sum=float(cfg["analysis"]["cpm_scale"]))
    sc.pp.log1p(adata)

    # Global perturbed vs control.
    sc.tl.rank_genes_groups(
        adata,
        groupby="de_global_group",
        groups=["perturbed"],
        reference="control",
        method=method,
        pts=True,
    )
    global_df = sc.get.rank_genes_groups_df(adata, group="perturbed")
    global_df.to_csv(out_scanpy / "de_global_wilcoxon.csv", index=False)
    safe_volcano(
        global_df,
        out_fig / "scanpy_global_wilcoxon_volcano.png",
        "Scanpy Wilcoxon: perturbed vs control",
        alpha,
    )

    # Per-target baseline.
    target_counts = (
        adata.obs.loc[adata.obs["parsed_target"] != "non-targeting", "parsed_target"]
        .astype(str)
        .value_counts()
    )
    eligible_targets = target_counts[target_counts >= min_cells_target].index.tolist()

    per_target_rows = []
    fail_rows = []
    per_target_out = out_scanpy / "de_per_target_wilcoxon.csv"
    existing = pd.DataFrame()
    done_targets = set()
    if per_target_out.exists():
        try:
            existing = pd.read_csv(per_target_out)
            if "grna_target" in existing.columns:
                done_targets = set(existing["grna_target"].astype(str).unique().tolist())
        except Exception:
            existing = pd.DataFrame()

    if len(existing) > 0:
        per_target_rows.append(existing)

    for tgt in eligible_targets:
        if tgt in done_targets:
            continue
        labels = np.where(
            adata.obs[control_col].astype(bool).values,
            "control",
            np.where(adata.obs["parsed_target"].astype(str).values == tgt, tgt, "other"),
        )
        if np.sum(labels == tgt) < min_cells_target:
            continue

        adata.obs["de_target_group"] = pd.Categorical(labels)
        try:
            sc.tl.rank_genes_groups(
                adata,
                groupby="de_target_group",
                groups=[tgt],
                reference="control",
                method=method,
                pts=True,
                n_genes=top_n_export,
            )
            df_t = sc.get.rank_genes_groups_df(adata, group=tgt)
            df_t["grna_target"] = tgt
            df_t = df_t.sort_values("pvals_adj", ascending=True).head(top_n_export)
            per_target_rows.append(df_t)
            pd.concat(per_target_rows, ignore_index=True).to_csv(per_target_out, index=False)
        except Exception as e:
            fail_rows.append({"grna_target": tgt, "error": str(e)})

    per_target_df = pd.concat(per_target_rows, ignore_index=True) if per_target_rows else pd.DataFrame()
    per_target_df.to_csv(per_target_out, index=False)

    if len(per_target_df) > 0:
        top_summary = (
            per_target_df.assign(significant=per_target_df["pvals_adj"] < alpha)
            .groupby("grna_target", as_index=False)
            .agg(
                n_genes=("names", "count"),
                n_sig=("significant", "sum"),
                best_gene=("names", "first"),
                best_padj=("pvals_adj", "min"),
            )
            .sort_values(["n_sig", "best_padj"], ascending=[False, True])
        )
    else:
        top_summary = pd.DataFrame(columns=["grna_target", "n_genes", "n_sig", "best_gene", "best_padj"])

    top_summary.to_csv(out_scanpy / "de_top_hits_summary.csv", index=False)
    pd.DataFrame(fail_rows).to_csv(out_scanpy / "de_per_target_failures.csv", index=False)


if __name__ == "__main__":
    main()
