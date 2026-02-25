# Pipeline Guide (for full-data runs)

This folder contains the runnable analysis pipeline. Canonical perturbation inference comes from **R `sceptre`**, then downstream analyses are run in Python.

## Pipeline stages

1. **R canonical inference**
- Script: `r_sceptre/run_sceptre_on_target_r.R`
- Produces on-target discovery table in `r_sceptre/outputs/`.

2. **Mirror R outputs to canonical location**
- Script: `run_all_r_canonical.sh` handles this automatically.
- Canonical tables go into `sceptre/`.

3. **Python downstream analyses**
- `run_on_target_expression_aware.py`
- `run_scanpy_wilcoxon_baseline.py`
- `run_scanpy_umap_leiden.py`
- `run_qc_marker_panels.py`
- `run_perturbation_burden.py`
- `run_pathway_hallmark.py`

## Files in this folder

- `config.yaml`: runtime config (paths, thresholds, clustering/DE settings)
- `run_all_r_canonical.sh`: one-command canonical run
- `run_all.sh`: legacy orchestration (not canonical source)
- stage scripts listed above
- `r_sceptre/run_sceptre_on_target_r.R`: R stage

## Required dependencies

### R
- `sceptre`
- `Matrix`

### Python
- `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `statsmodels`, `gseapy`, `pyyaml`

## Before you run

Edit `config.yaml`:
- set `data.h5ad_path` to your full dataset path,
- set `data.results_root` to your desired output location,
- set `pathways.hallmark_gmt` to a local Hallmark GMT file.

## Run

```bash
bash run_all_r_canonical.sh
```

## Outputs

- `sceptre/`: canonical R-based inference + expression-aware filtered tables
- `scanpy/`: DE baseline, embeddings, clustering summaries
- `pathways/`: Hallmark enrichment outputs and run log
- `perturbation_burden/`: per-cell burden tables
- `figures/`: plots

## Notes

- Controls are defined as `control_cell == True`.
- Guide assignment is derived from `assigned_guide` labels for this dataset format.
- For full-scale data, run on server-local storage (not Samba-mounted compute paths).
