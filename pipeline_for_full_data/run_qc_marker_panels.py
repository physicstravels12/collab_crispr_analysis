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


def res_tag(r: float) -> str:
    return str(r).replace(".", "p")


def scatter_cont(df: pd.DataFrame, value_col: str, out_png: Path, title: str):
    if value_col not in df.columns:
        return
    plt.figure(figsize=(7, 6))
    v = pd.to_numeric(df[value_col], errors="coerce")
    idx = np.isfinite(v.values)
    sca = plt.scatter(df.loc[idx, "UMAP1"], df.loc[idx, "UMAP2"], c=v[idx], s=3, alpha=0.8, cmap="viridis", linewidths=0)
    plt.colorbar(sca, label=value_col)
    plt.title(title)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def main():
    root = Path(__file__).resolve().parent
    cfg = yaml.safe_load((root / "config.yaml").read_text(encoding="utf-8"))
    if not bool(cfg.get("extras_qc_markers", {}).get("enabled", True)):
        return

    out_scanpy = root / "scanpy"
    out_fig = root / "figures"
    out_scanpy.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)

    emb_path = out_scanpy / "embeddings_umap.csv"
    if not emb_path.exists():
        raise FileNotFoundError(f"Missing {emb_path}. Run run_scanpy_umap_leiden.py first.")

    emb = pd.read_csv(emb_path)
    adata = ad.read_h5ad(cfg["data"]["h5ad_path"])
    adata = adata.copy()

    cell_to_cluster = emb.set_index("cell_id")
    default_res = float(cfg["scanpy_cluster"]["default_resolution"])
    cluster_col = f"leiden_{res_tag(default_res)}"
    if cluster_col not in cell_to_cluster.columns:
        raise RuntimeError(f"Default cluster column missing: {cluster_col}")

    cluster_vals = cell_to_cluster.loc[adata.obs_names.astype(str), cluster_col].astype(str).values
    adata.obs[cluster_col] = pd.Categorical(cluster_vals)

    # Marker genes per cluster.
    sc.pp.normalize_total(adata, target_sum=float(cfg["analysis"]["cpm_scale"]))
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, groupby=cluster_col, method="wilcoxon", pts=True)

    top_n = int(cfg["extras_qc_markers"].get("top_markers_per_cluster", 20))
    marker_rows = []
    for c in sorted(adata.obs[cluster_col].astype(str).unique().tolist()):
        df = sc.get.rank_genes_groups_df(adata, group=c).head(top_n)
        df["cluster"] = c
        marker_rows.append(df)
    markers = pd.concat(marker_rows, ignore_index=True) if marker_rows else pd.DataFrame()
    markers.to_csv(out_scanpy / "cluster_markers_topN.csv", index=False)

    # Composition control vs perturbed by cluster.
    control_col = cfg["analysis"]["control_column"]
    comp = (
        pd.DataFrame({
            "cluster": adata.obs[cluster_col].astype(str).values,
            "control_cell": adata.obs[control_col].astype(bool).values,
        })
        .groupby(["cluster", "control_cell"], as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
    )
    comp["frac_within_cluster"] = comp["n_cells"] / comp.groupby("cluster")["n_cells"].transform("sum")
    comp.to_csv(out_scanpy / "cluster_composition_control_vs_perturbed.csv", index=False)

    # Composition by parsed target (top targets + other).
    if "parsed_target" in emb.columns:
        min_n = int(cfg["extras_qc_markers"].get("min_cells_per_cluster_composition", 20))
        tmp = emb[["cell_id", "parsed_target", cluster_col]].copy()
        target_counts = tmp.loc[tmp["parsed_target"] != "non-targeting", "parsed_target"].value_counts()
        top_targets = set(target_counts.head(12).index.tolist())
        tmp["target_group"] = np.where(
            tmp["parsed_target"] == "non-targeting",
            "control_nt",
            np.where(tmp["parsed_target"].isin(top_targets), tmp["parsed_target"], "other_perturbed"),
        )
        comp_t = (
            tmp.groupby([cluster_col, "target_group"], as_index=False)
            .size()
            .rename(columns={"size": "n_cells", cluster_col: "cluster"})
        )
        comp_t = comp_t[comp_t["n_cells"] >= min_n]
        comp_t["frac_within_cluster"] = comp_t["n_cells"] / comp_t.groupby("cluster")["n_cells"].transform("sum")
        comp_t.to_csv(out_scanpy / "cluster_composition_by_target.csv", index=False)

    # QC overlays.
    qcdf = emb.copy()
    burden_path = root / "perturbation_burden" / "per_cell_scores.csv"
    if burden_path.exists():
        burden = pd.read_csv(burden_path)
        if "perturbation_burden_weighted" in burden.columns and "cell_id" in burden.columns:
            qcdf = qcdf.merge(
                burden[["cell_id", "perturbation_burden_weighted"]],
                on="cell_id",
                how="left",
            )

    scatter_cont(qcdf, "num_umis", out_fig / "umap_qc_num_umis.png", "UMAP QC: num_umis")
    scatter_cont(qcdf, "num_features", out_fig / "umap_qc_num_features.png", "UMAP QC: num_features")
    scatter_cont(
        qcdf,
        "perturbation_burden_weighted",
        out_fig / "umap_qc_burden.png",
        "UMAP QC: perturbation burden",
    )


if __name__ == "__main__":
    main()
