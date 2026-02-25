# CRISPR Perturb-seq Analysis (R-canonical)

This repository packages two things in one place:

- `results_snapshot_current_data/`: a frozen snapshot of analysis outputs from the current subset dataset.
- `pipeline_for_full_data/`: the runnable pipeline collaborators should use on larger/full datasets.

The canonical perturbation inference source in this repo is **R `sceptre`** (not the deprecated Python port).

## Why this repo exists
This project analyzes single-cell CRISPR perturbation data with emphasis on:
- on-target perturbation effects,
- global/per-target DE baselines,
- UMAP/Leiden cell-state structure,
- pathway-level shifts (Hallmark),
- per-cell perturbation burden.

A previous Python `sceptre` route showed directionality inconsistencies for on-target fold-changes in this data context. For that reason, this repo standardizes on **R `sceptre`** output, then performs downstream analyses in Python.

## Repository structure

```text
.
├── pipeline_for_full_data/
│   ├── config.yaml
│   ├── run_all_r_canonical.sh
│   ├── run_scanpy_wilcoxon_baseline.py
│   ├── run_scanpy_umap_leiden.py
│   ├── run_qc_marker_panels.py
│   ├── run_on_target_expression_aware.py
│   ├── run_perturbation_burden.py
│   ├── run_pathway_hallmark.py
│   ├── r_sceptre/run_sceptre_on_target_r.R
│   └── README.md
└── results_snapshot_current_data/
    ├── analysis_summary.txt
    ├── sceptre/
    ├── scanpy/
    ├── pathways/
    ├── perturbation_burden/
    └── figures/
```

## Upstream tools and repos this relies on

1. **R `sceptre` package** (Katsevich Lab)
- Upstream repo: https://github.com/Katsevich-Lab/sceptre
- Used for canonical perturbation discovery statistics.

2. **Python scientific stack**
- `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `statsmodels`, `gseapy`, `pyyaml`
- Used for downstream DE baselines, clustering, pathway, and burden analyses.

3. **MSigDB Hallmark gene sets**
- Pipeline expects a local Hallmark GMT path in config (`pathways.hallmark_gmt`).

## Data assumptions / input contract

The pipeline expects an `.h5ad` with:
- expression matrix in `adata.X` (cells x genes),
- `obs` columns:
  - `control_cell` (boolean),
  - `assigned_guide` (string per cell),
- gene names in `adata.var_names`.

Guide names are parsed into target labels by convention (`NO_SITE` / `NON-GENE_SITE` -> `non-targeting`).

## Canonical workflow (recommended)

Run from `pipeline_for_full_data/`:

```bash
bash run_all_r_canonical.sh
```

This does:
1. Run R `sceptre` (`r_sceptre/run_sceptre_on_target_r.R`).
2. Mirror R outputs into `sceptre/` canonical CSVs.
3. Run Python downstream analyses:
   - expression-aware on-target filtering,
   - Scanpy Wilcoxon baseline,
   - UMAP + Leiden,
   - QC/marker panels,
   - perturbation burden,
   - Hallmark pathway enrichment.

## What to configure before running on full data

Edit `pipeline_for_full_data/config.yaml`:
- `data.h5ad_path`: absolute server path to full `.h5ad`.
- `data.results_root`: desired output root on server storage.
- `pathways.hallmark_gmt`: absolute path to local Hallmark GMT file.

Important: the checked-in config currently contains local development paths and must be updated for your environment.

## Output map

Main outputs are written under the configured `results_root`:

- `sceptre/`
  - canonical R-based on-target/discovery tables
  - expression-aware filtered tables
- `scanpy/`
  - global/per-target DE baseline tables
  - UMAP embeddings and cluster-level summary tables
- `pathways/`
  - Hallmark enrichment outputs/log
- `perturbation_burden/`
  - per-cell perturbation burden tables
- `figures/`
  - all plots

## Recommended environment / execution model

For full datasets (e.g., 250GB), run on institute server close to data:
- connect via SSH (optionally VS Code Remote SSH),
- run in persistent shell (`tmux`/`nohup`),
- avoid compute over Samba mounts.

## Validation checks after a run

1. Canonical source confirmation:
- `sceptre/result_source.csv` should indicate `R_sceptre`.

2. Expected key files exist:
- `sceptre/on_target_results.csv`
- `scanpy/de_global_wilcoxon.csv`
- `scanpy/de_per_target_wilcoxon.csv`
- `scanpy/embeddings_umap.csv`
- `pathways/hallmark_global.csv`
- `perturbation_burden/per_cell_scores.csv`

3. Review `analysis_summary.txt` for quick run-level interpretation.

## Snapshot vs pipeline

- `results_snapshot_current_data/` is historical output and should not be edited for new runs.
- `pipeline_for_full_data/` is the operational code to run on new/full datasets.

## Known caveats

- Hallmark enrichment requires a valid local GMT file path.
- Target-level Hallmark can be empty if each target contributes too few ranked genes (e.g., on-target-only single-gene test sets).
- If running outside development environment, update path assumptions in config before execution.
