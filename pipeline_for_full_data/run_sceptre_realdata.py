#!/usr/bin/env python3
import json
import os
import time
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

import sceptre


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


def mean_var_sparse(X: sp.spmatrix):
    mean = np.asarray(X.mean(axis=0)).ravel()
    mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
    var = np.maximum(mean_sq - mean ** 2, 0)
    return mean, var


def log1p_cpm_sparse(X_csr: sp.csr_matrix, cpm_scale: float) -> sp.csr_matrix:
    libsize = np.asarray(X_csr.sum(axis=1)).ravel().astype(np.float64)
    scale = np.where(libsize > 0, cpm_scale / libsize, 0.0)
    X_cpm = X_csr.multiply(scale[:, None]).tocsr().astype(np.float32)
    X_log = X_cpm.copy()
    X_log.data = np.log1p(X_log.data)
    return X_log


def run_batch(
    response_matrix,
    grna_matrix,
    grna_target_df,
    response_names,
    grna_names,
    extra_covariates,
    pairs_df,
    batch_name,
    out_csv,
    cfg,
):
    obj = sceptre.import_data(
        response_matrix=response_matrix,
        grna_matrix=grna_matrix,
        grna_target_data_frame=grna_target_df,
        moi=cfg["sceptre"]["moi"],
        extra_covariates=extra_covariates,
        response_names=response_names,
        grna_names=grna_names,
    )

    obj.set_analysis_parameters(
        discovery_pairs=pairs_df,
        side=cfg["sceptre"]["side"],
        control_group=cfg["sceptre"]["control_group"],
        resampling_mechanism=cfg["sceptre"]["resampling_mechanism"],
        multiple_testing_method=cfg["sceptre"]["multiple_testing_method"],
        multiple_testing_alpha=cfg["sceptre"]["multiple_testing_alpha"],
    )

    assign_method = cfg["sceptre"]["assign_method"]
    if assign_method == "thresholding":
        obj.assign_grnas(method=assign_method, threshold=float(cfg["sceptre"]["assign_threshold"]))
    else:
        obj.assign_grnas(method=assign_method)
    obj.run_qc(**cfg["sceptre"]["run_qc"])

    obj.B1 = int(cfg["sceptre"]["B1"])
    obj.B2 = int(cfg["sceptre"]["B2"])

    obj.run_discovery_analysis(print_progress=False)
    res = obj.get_result("discovery").copy()
    res["batch_name"] = batch_name
    res.to_csv(out_csv, index=False)
    return {
        "batch_name": batch_name,
        "n_pairs_input": int(len(pairs_df)),
        "n_pairs_output": int(len(res)),
        "n_pass_qc": int(res["pass_qc"].sum()) if "pass_qc" in res.columns else None,
    }


def main():
    cfg_path = Path(__file__).resolve().parent / "config.yaml"
    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    np.random.seed(int(cfg["analysis"]["seed"]))

    results_root = Path(cfg["data"]["results_root"])
    out_sceptre = results_root / "sceptre"
    out_fig = results_root / "figures"
    out_sceptre.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(cfg["data"]["h5ad_path"])
    obs = adata.obs.copy()
    var_names = adata.var_names.astype(str).tolist()

    control_col = cfg["analysis"]["control_column"]
    guide_col = cfg["analysis"]["assigned_guide_column"]

    obs[guide_col] = obs[guide_col].astype(str)
    obs["parsed_target"] = obs[guide_col].map(parse_target_from_guide)

    ctrl_mask = obs[control_col].astype(bool).values
    ctrl_targets = set(obs.loc[ctrl_mask, "parsed_target"].astype(str).unique())
    if ctrl_targets != {"non-targeting"}:
        raise RuntimeError(f"Expected control cells to map only to non-targeting, saw: {sorted(ctrl_targets)}")

    guides = sorted(obs[guide_col].astype(str).unique().tolist())
    guide_to_idx = {g: i for i, g in enumerate(guides)}
    guide_idx = obs[guide_col].map(guide_to_idx).values

    n_cells = adata.n_obs
    n_guides = len(guides)
    data = np.ones(n_cells, dtype=np.float32)
    cols = np.arange(n_cells, dtype=np.int64)
    grna_matrix = sp.csr_matrix((data, (guide_idx, cols)), shape=(n_guides, n_cells), dtype=np.float32)

    sums = np.asarray(grna_matrix.sum(axis=0)).ravel()
    if not np.all(sums == 1):
        raise RuntimeError("Derived gRNA matrix is not one-hot per cell")

    guide_targets = [parse_target_from_guide(g) for g in guides]
    grna_target_df = pd.DataFrame(
        {
            "grna_id": guides,
            "grna_target": guide_targets,
            "chr": [""] * len(guides),
            "start": [0] * len(guides),
            "end": [0] * len(guides),
        }
    )

    if sp.issparse(adata.X):
        X_cells_genes = adata.X.tocsr()
    else:
        X_cells_genes = sp.csr_matrix(adata.X)
    response_matrix = X_cells_genes.T.tocsr().astype(np.float32)

    extra_covariates = pd.DataFrame(index=np.arange(n_cells))

    target_counts = (
        obs.loc[obs["parsed_target"] != "non-targeting", "parsed_target"]
        .astype(str)
        .value_counts()
    )
    eligible_targets = target_counts[target_counts >= int(cfg["analysis"]["min_cells_per_target"])].index.tolist()

    X_log = log1p_cpm_sparse(X_cells_genes, float(cfg["analysis"]["cpm_scale"]))
    detect = np.asarray((X_cells_genes > 0).mean(axis=0)).ravel()
    _, var_log = mean_var_sparse(X_log)

    keep_genes = detect >= float(cfg["analysis"]["min_detection_fraction"])
    keep_gene_set = {var_names[i] for i in np.where(keep_genes)[0]}
    panel_size = int(cfg["analysis"]["trans_panel_size"])
    rank_idx = np.argsort(-var_log)
    trans_genes = [var_names[i] for i in rank_idx if keep_genes[i]][:panel_size]

    discovery_mode = str(cfg["analysis"].get("discovery_mode", "combined"))

    on_target_pairs = []
    for tgt in eligible_targets:
        if tgt in set(var_names) and tgt in keep_gene_set:
            on_target_pairs.append({"grna_target": tgt, "response_id": tgt, "pair_type": "on_target"})

    trans_pairs = []
    if discovery_mode != "on_target_only":
        for tgt in eligible_targets:
            for g in trans_genes:
                trans_pairs.append({"grna_target": tgt, "response_id": g, "pair_type": "trans_panel"})

    pairs_all = pd.DataFrame(trans_pairs + on_target_pairs)
    pairs_all = pairs_all.drop_duplicates(subset=["grna_target", "response_id"]).reset_index(drop=True)

    pairs_all.to_csv(out_sceptre / "discovery_pairs_all.csv", index=False)
    pd.Series(eligible_targets, name="target").to_csv(out_sceptre / "eligible_targets.csv", index=False)
    pd.Series(trans_genes, name="gene").to_csv(out_sceptre / "trans_panel_genes.csv", index=False)

    pilot_n = int(cfg["execution"]["pilot_n_targets"])
    targets_sorted = target_counts.loc[eligible_targets].sort_values(ascending=False).index.tolist()
    pilot_targets = targets_sorted[:pilot_n]
    full_targets = targets_sorted[pilot_n:]

    run_metrics = {
        "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_cells": int(n_cells),
        "n_genes": int(len(var_names)),
        "n_guides": int(n_guides),
        "n_eligible_targets": int(len(eligible_targets)),
        "n_on_target_pairs": int(len(on_target_pairs)),
        "n_trans_pairs": int(len(trans_pairs)),
        "n_total_pairs": int(len(pairs_all)),
        "batches": [],
        "errors": [],
    }

    single_run = bool(cfg["execution"].get("single_run", False))
    if single_run:
        one_csv = out_sceptre / "batch_single.csv"
        t0 = time.time()
        try:
            if not one_csv.exists():
                m = run_batch(
                    response_matrix,
                    grna_matrix,
                    grna_target_df,
                    var_names,
                    guides,
                    extra_covariates,
                    pairs_all[["grna_target", "response_id"]].copy(),
                    "single",
                    one_csv,
                    cfg,
                )
            else:
                m = {"batch_name": "single", "skipped_existing": True, "n_pairs_input": int(len(pairs_all))}
            m["runtime_sec"] = round(time.time() - t0, 2)
            run_metrics["batches"].append(m)
        except Exception as e:
            run_metrics["errors"].append({"batch": "single", "error": str(e)})
    else:
        pilot_pairs = pairs_all[pairs_all["grna_target"].isin(pilot_targets)][["grna_target", "response_id"]].copy()
        pilot_out = out_sceptre / "batch_pilot.csv"
        t0 = time.time()
        try:
            if not pilot_out.exists():
                m = run_batch(
                    response_matrix,
                    grna_matrix,
                    grna_target_df,
                    var_names,
                    guides,
                    extra_covariates,
                    pilot_pairs,
                    "pilot",
                    pilot_out,
                    cfg,
                )
            else:
                m = {"batch_name": "pilot", "skipped_existing": True, "n_pairs_input": int(len(pilot_pairs))}
            m["runtime_sec"] = round(time.time() - t0, 2)
            run_metrics["batches"].append(m)
        except Exception as e:
            run_metrics["errors"].append({"batch": "pilot", "error": str(e)})

        if bool(cfg["execution"].get("run_full", True)):
            bs = int(cfg["execution"]["batch_size_targets"])
            for bi, start in enumerate(range(0, len(full_targets), bs), 1):
                target_chunk = full_targets[start : start + bs]
                batch_name = f"full_{bi:03d}"
                out_csv = out_sceptre / f"batch_{batch_name}.csv"
                batch_pairs = pairs_all[pairs_all["grna_target"].isin(target_chunk)][["grna_target", "response_id"]].copy()
                t1 = time.time()
                try:
                    if out_csv.exists():
                        m = {
                            "batch_name": batch_name,
                            "skipped_existing": True,
                            "n_pairs_input": int(len(batch_pairs)),
                            "n_targets": int(len(target_chunk)),
                        }
                    else:
                        m = run_batch(
                            response_matrix,
                            grna_matrix,
                            grna_target_df,
                            var_names,
                            guides,
                            extra_covariates,
                            batch_pairs,
                            batch_name,
                            out_csv,
                            cfg,
                        )
                        m["n_targets"] = int(len(target_chunk))
                    m["runtime_sec"] = round(time.time() - t1, 2)
                    run_metrics["batches"].append(m)
                except Exception as e:
                    run_metrics["errors"].append({"batch": batch_name, "error": str(e)})

    batch_files = sorted(out_sceptre.glob("batch_*.csv"))
    if batch_files:
        merged = pd.concat([pd.read_csv(p) for p in batch_files], ignore_index=True)
        merged = merged.drop_duplicates(subset=["grna_target", "response_id"], keep="first")
        merged.to_csv(out_sceptre / "discovery_merged.csv", index=False)

        on_target = merged[merged["grna_target"].astype(str) == merged["response_id"].astype(str)].copy()
        trans_only = merged[merged["grna_target"].astype(str) != merged["response_id"].astype(str)].copy()

        on_target.to_csv(out_sceptre / "on_target_results.csv", index=False)
        trans_only.to_csv(out_sceptre / "trans_panel_results.csv", index=False)

        if len(on_target) > 0 and "log_2_fold_change" in on_target.columns:
            s = on_target.sort_values("log_2_fold_change", ascending=True).head(30)
            plt.figure(figsize=(9, 5))
            plt.barh(s["grna_target"], s["log_2_fold_change"], color="#2563eb")
            plt.axvline(0, color="black", linewidth=1)
            plt.xlabel("on-target log2 fold-change")
            plt.ylabel("target")
            plt.title("Top 30 strongest negative on-target effects")
            plt.tight_layout()
            plt.savefig(out_fig / "on_target_top30_negative_log2fc.png", dpi=180)
            plt.close()

        if "p_value" in merged.columns:
            p = np.clip(merged["p_value"].astype(float).values, 1e-300, 1.0)
            plt.figure(figsize=(6, 4))
            plt.hist(p, bins=40, color="#374151", alpha=0.9)
            plt.xlabel("p-value")
            plt.ylabel("count")
            plt.title("SCEPTRE discovery p-value distribution")
            plt.tight_layout()
            plt.savefig(out_fig / "sceptre_pvalue_histogram.png", dpi=180)
            plt.close()

        run_metrics["n_merged_results"] = int(len(merged))
        run_metrics["n_on_target_results"] = int(len(on_target))

    run_metrics["finished_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(out_sceptre / "pilot_metrics.json", "w", encoding="utf-8") as f:
        json.dump(run_metrics, f, indent=2)

    summary_lines = [
        "SCEPTRE real-data run summary",
        f"cells={run_metrics['n_cells']}, genes={run_metrics['n_genes']}, guides={run_metrics['n_guides']}",
        f"eligible_targets={run_metrics['n_eligible_targets']}",
        f"total_pairs={run_metrics['n_total_pairs']} (trans={run_metrics['n_trans_pairs']}, on_target={run_metrics['n_on_target_pairs']})",
        f"batches_executed_or_skipped={len(run_metrics['batches'])}",
        f"errors={len(run_metrics['errors'])}",
    ]
    if run_metrics.get("n_merged_results") is not None:
        summary_lines.append(f"merged_results={run_metrics['n_merged_results']}")
    if run_metrics.get("n_on_target_results") is not None:
        summary_lines.append(f"on_target_results={run_metrics['n_on_target_results']}")

    (out_sceptre / "run_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
