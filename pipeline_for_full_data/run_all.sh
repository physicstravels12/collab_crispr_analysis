#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/ninaanikeeva/Desktop/MIT/research/bio/results1/sceptre_realdata"

python3 "$ROOT/run_sceptre_realdata.py"
python3 "$ROOT/run_scanpy_wilcoxon_baseline.py"
python3 "$ROOT/run_scanpy_umap_leiden.py"
python3 "$ROOT/run_on_target_expression_aware.py"
python3 "$ROOT/run_qc_marker_panels.py"
python3 "$ROOT/run_pathway_hallmark.py"
python3 "$ROOT/run_perturbation_burden.py"

echo "Completed: $ROOT"
