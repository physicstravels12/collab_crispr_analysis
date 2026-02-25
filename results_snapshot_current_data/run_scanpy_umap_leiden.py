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


def res_tag(r: float) -> str:
    return str(r).replace(".", "p")


def plot_umap_categorical(df: pd.DataFrame, color_col: str, title: str, out_png: Path, max_categories: int = 25):
    vals = df[color_col].astype(str)
    cats = vals.value_counts().index.tolist()
    if len(cats) > max_categories:
        keep = set(cats[: max_categories - 1])
        vals = vals.where(vals.isin(keep), "other")
        cats = vals.value_counts().index.tolist()
    cmap = plt.get_cmap("tab20")
    colors = {c: cmap(i % 20) for i, c in enumerate(cats)}

    plt.figure(figsize=(7, 6))
    for c in cats:
        m = vals == c
        plt.scatter(df.loc[m, "UMAP1"], df.loc[m, "UMAP2"], s=3, alpha=0.7, c=[colors[c]], label=c, linewidths=0)
    plt.title(title)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.legend(markerscale=3, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def main():
    root = Path(__file__).resolve().parent
    cfg = yaml.safe_load((root / "config.yaml").read_text(encoding="utf-8"))
    if not bool(cfg.get("scanpy_cluster", {}).get("enabled", True)):
        return

    out_scanpy = root / "scanpy"
    out_fig = root / "figures"
    out_scanpy.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(cfg["data"]["h5ad_path"])
    adata = adata.copy()

    control_col = cfg["analysis"]["control_column"]
    guide_col = cfg["analysis"]["assigned_guide_column"]
    seed = int(cfg["analysis"]["seed"])

    adata.obs[guide_col] = adata.obs[guide_col].astype(str)
    adata.obs["parsed_target"] = adata.obs[guide_col].map(parse_target_from_guide)
    adata.obs["perturbation_group"] = np.where(adata.obs[control_col].astype(bool), "control", "perturbed")

    sc.pp.normalize_total(adata, target_sum=float(cfg["analysis"]["cpm_scale"]))
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=int(cfg["scanpy_cluster"]["n_hvgs"]), subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=int(cfg["scanpy_cluster"]["n_pcs"]), random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=int(cfg["scanpy_cluster"]["neighbors_k"]), n_pcs=int(cfg["scanpy_cluster"]["n_pcs"]))
    sc.tl.umap(
        adata,
        min_dist=float(cfg["scanpy_cluster"]["umap_min_dist"]),
        spread=float(cfg["scanpy_cluster"]["umap_spread"]),
        random_state=seed,
    )

    resolutions = [float(x) for x in cfg["scanpy_cluster"]["leiden_resolutions"]]
    leiden_cols = []
    size_rows = []
    for r in resolutions:
        key = f"leiden_{res_tag(r)}"
        sc.tl.leiden(adata, resolution=r, key_added=key, random_state=seed)
        leiden_cols.append(key)
        counts = adata.obs[key].value_counts().sort_index()
        for cluster, n in counts.items():
            size_rows.append({"resolution": r, "cluster": str(cluster), "n_cells": int(n)})

    umap = adata.obsm["X_umap"]
    emb = pd.DataFrame({"cell_id": adata.obs_names.astype(str), "UMAP1": umap[:, 0], "UMAP2": umap[:, 1]})
    keep_obs = [control_col, guide_col, "parsed_target", "perturbation_group", "num_umis", "num_features", "lane", "sample_id"]
    keep_obs = [c for c in keep_obs if c in adata.obs.columns]
    emb = emb.join(adata.obs[keep_obs].reset_index(drop=True))
    for c in leiden_cols:
        emb[c] = adata.obs[c].astype(str).values

    emb.to_csv(out_scanpy / "embeddings_umap.csv", index=False)
    pd.DataFrame(size_rows).to_csv(out_scanpy / "leiden_cluster_sizes.csv", index=False)

    for r in resolutions:
        c = f"leiden_{res_tag(r)}"
        plot_umap_categorical(
            emb,
            c,
            f"UMAP Leiden clusters (resolution={r})",
            out_fig / f"umap_leiden_res_{res_tag(r)}.png",
            max_categories=40,
        )

    plot_umap_categorical(
        emb,
        "perturbation_group",
        "UMAP: control vs perturbed",
        out_fig / "umap_control_vs_perturbed.png",
        max_categories=5,
    )

    top_n_targets = int(cfg["analysis"].get("top_targets_umap_label", 10))
    target_counts = emb.loc[emb["parsed_target"] != "non-targeting", "parsed_target"].value_counts()
    top_targets = set(target_counts.head(top_n_targets).index.tolist())
    tgt_label = emb["parsed_target"].astype(str)
    tgt_label = np.where(tgt_label == "non-targeting", "control_nt", np.where(pd.Series(tgt_label).isin(top_targets), tgt_label, "other_perturbed"))
    emb["target_top_label"] = tgt_label
    plot_umap_categorical(
        emb,
        "target_top_label",
        f"UMAP: top {top_n_targets} targets + control",
        out_fig / "umap_by_top_targets.png",
        max_categories=20,
    )


if __name__ == "__main__":
    main()
