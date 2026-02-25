#!/usr/bin/env python3
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/.cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

import gseapy as gp


def run_prerank_safe(rnk: pd.DataFrame, gene_sets, seed: int):
    return gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=None,
        min_size=10,
        max_size=500,
        permutation_num=500,
        seed=seed,
        verbose=False,
    )


def main():
    root = Path(__file__).resolve().parent
    cfg = yaml.safe_load((root / "config.yaml").read_text(encoding="utf-8"))

    out_dir = root / "pathways"
    fig_dir = root / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    disc_path = root / "sceptre" / "discovery_merged.csv"
    if not disc_path.exists():
        raise FileNotFoundError(f"Missing {disc_path}")

    disc = pd.read_csv(disc_path)
    if len(disc) == 0:
        (out_dir / "hallmark_target_level.csv").write_text("", encoding="utf-8")
        (out_dir / "hallmark_global.csv").write_text("", encoding="utf-8")
        return

    gene_sets = cfg["pathways"].get("hallmark_gmt")
    if gene_sets is not None:
        gene_sets = str(gene_sets)

    log = {"gene_sets": gene_sets, "errors": []}
    if gene_sets is None:
        gene_sets = "Hallmark_2020"

    target_results = []
    seed = int(cfg["analysis"]["seed"])

    # Per-target Hallmark using SCEPTRE-ranked genes.
    try:
        for tgt, sub in disc.groupby("grna_target"):
            tmp = sub.copy()
            p = np.clip(tmp["p_value"].astype(float).values, 1e-300, 1.0)
            sign = np.sign(tmp["log_2_fold_change"].astype(float).values)
            score = sign * -np.log10(p)
            rnk = pd.DataFrame({"gene": tmp["response_id"].astype(str).values, "score": score})
            rnk = rnk.drop_duplicates(subset=["gene"], keep="first")
            if len(rnk) < 30:
                continue
            pre = run_prerank_safe(rnk, gene_sets, seed)
            res = pre.res2d.reset_index().rename(columns={"Term": "pathway"})
            if len(res) == 0:
                continue
            res["grna_target"] = tgt
            target_results.append(res)
    except Exception as e:
        log["errors"].append(f"target_level_failed: {e}")

    if target_results:
        target_df = pd.concat(target_results, ignore_index=True)
        target_df.to_csv(out_dir / "hallmark_target_level.csv", index=False)

        # Heatmap of top pathways by mean NES magnitude.
        top_n = int(cfg["pathways"]["top_n_pathways_plot"])
        agg = (
            target_df.groupby("pathway", as_index=False)["NES"]
            .apply(lambda s: np.mean(np.abs(s)))
            .rename(columns={"NES": "mean_abs_nes"})
            .sort_values("mean_abs_nes", ascending=False)
            .head(top_n)
        )
        keep_paths = agg["pathway"].tolist()
        mat = (
            target_df[target_df["pathway"].isin(keep_paths)]
            .pivot_table(index="pathway", columns="grna_target", values="NES", aggfunc="first")
            .fillna(0.0)
        )
        if mat.shape[0] > 0 and mat.shape[1] > 0:
            plt.figure(figsize=(max(8, 0.2 * mat.shape[1]), max(5, 0.35 * mat.shape[0])))
            sns.heatmap(mat, cmap="coolwarm", center=0, cbar_kws={"label": "NES"})
            plt.title("Hallmark pathway NES across perturbation targets")
            plt.xlabel("target")
            plt.ylabel("pathway")
            plt.tight_layout()
            plt.savefig(fig_dir / "hallmark_target_heatmap.png", dpi=180)
            plt.close()
    else:
        pd.DataFrame().to_csv(out_dir / "hallmark_target_level.csv", index=False)

    # Global Hallmark from existing global DE table.
    global_de = Path("/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/de_global_perturbed_vs_controls.csv")
    try:
        g = pd.read_csv(global_de)
        p = np.clip(g["p_value"].astype(float).values, 1e-300, 1.0)
        sign = np.sign(g["log2fc"].astype(float).values)
        score = sign * -np.log10(p)
        rnk = pd.DataFrame({"gene": g["gene"].astype(str).values, "score": score})
        rnk = rnk.drop_duplicates(subset=["gene"], keep="first")
        pre = run_prerank_safe(rnk, gene_sets, seed)
        glob = pre.res2d.reset_index().rename(columns={"Term": "pathway"})
        glob.to_csv(out_dir / "hallmark_global.csv", index=False)

        if len(glob) > 0:
            top = glob.sort_values("FDR q-val", ascending=True).head(15)
            plt.figure(figsize=(7, max(4, 0.3 * len(top))))
            plt.barh(top["pathway"], top["NES"], color="#0f766e")
            plt.axvline(0, color="black", linewidth=1)
            plt.xlabel("NES")
            plt.ylabel("pathway")
            plt.title("Global Hallmark enrichment (perturbed vs control)")
            plt.tight_layout()
            plt.savefig(fig_dir / "hallmark_global_top15.png", dpi=180)
            plt.close()
    except Exception as e:
        log["errors"].append(f"global_failed: {e}")
        pd.DataFrame().to_csv(out_dir / "hallmark_global.csv", index=False)

    (out_dir / "hallmark_run_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
