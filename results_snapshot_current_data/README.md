# Real-data SCEPTRE pipeline (`singlets.h5ad`)

This folder contains a reproducible analysis workflow for:
- on-target knockdown with R `sceptre` as canonical source,
- Scanpy Wilcoxon DE baseline (global + per-target),
- UMAP + Leiden clustering (resolution sweep),
- expression-aware on-target filtering,
- cluster markers and composition diagnostics,
- pathway shifts (Hallmark GSEA),
- per-cell perturbation burden scalar.

## Files
- `config.yaml`: all thresholds, run sizing, and output paths.
- `run_sceptre_realdata.py`: legacy non-canonical SCEPTRE discovery script.
- `run_scanpy_wilcoxon_baseline.py`: global and per-target Wilcoxon DE baseline.
- `run_scanpy_umap_leiden.py`: UMAP embedding + Leiden clustering.
- `run_on_target_expression_aware.py`: expression-aware on-target post-filtering.
- `run_qc_marker_panels.py`: cluster marker calling and composition/QC overlays.
- `run_pathway_hallmark.py`: Hallmark enrichment from SCEPTRE and global DE rankings.
- `run_perturbation_burden.py`: per-cell perturbation burden score.
- `run_all.sh`: orchestrates full pipeline.
- `run_all_r_canonical.sh`: orchestrates full pipeline with R sceptre as the canonical source.

## Outputs
- `sceptre/discovery_merged.csv`
- `sceptre/on_target_results.csv`
- `sceptre/on_target_results_expressed_only.csv`
- `sceptre/on_target_expression_filter_stats.csv`
- `sceptre/trans_panel_results.csv`
- `sceptre/pilot_metrics.json`
- `scanpy/de_global_wilcoxon.csv`
- `scanpy/de_per_target_wilcoxon.csv`
- `scanpy/embeddings_umap.csv`
- `scanpy/leiden_cluster_sizes.csv`
- `scanpy/cluster_markers_topN.csv`
- `scanpy/cluster_composition_control_vs_perturbed.csv`
- `scanpy/cluster_composition_by_target.csv`
- `pathways/hallmark_target_level.csv`
- `pathways/hallmark_global.csv`
- `perturbation_burden/per_cell_scores.csv`
- `figures/*.png`

## Run
```bash
bash /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_all.sh
```

Or stepwise:
```bash
python3 /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_sceptre_realdata.py
python3 /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_scanpy_wilcoxon_baseline.py
python3 /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_scanpy_umap_leiden.py
python3 /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_on_target_expression_aware.py
python3 /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_qc_marker_panels.py
python3 /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_pathway_hallmark.py
python3 /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_perturbation_burden.py
```

R-canonical one-command run:
```bash
bash /Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata/run_all_r_canonical.sh
```

## Notes
- Controls are defined as `control_cell == True` and are validated to map to parsed `non-targeting` guide families.
- Because `singlets.h5ad` stores assigned guide labels (not raw gRNA counts), the gRNA matrix is derived as one-hot guide-by-cell assignments.
- The current config uses `thresholding` assignment (`threshold=1`) to keep one-hot assignments consistent with low-MOI singlet data.
- Batches are checkpointed as `sceptre/batch_*.csv`; reruns skip completed batches automatically.
