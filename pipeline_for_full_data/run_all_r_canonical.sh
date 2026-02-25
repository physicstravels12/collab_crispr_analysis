#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$ROOT/../.." && pwd)"
RSCRIPT="${RSCRIPT_BIN:-$(command -v Rscript || true)}"
if [[ -z "$RSCRIPT" ]]; then
  echo "Rscript not found in PATH. Set RSCRIPT_BIN or install R."
  exit 1
fi
export R_LIBS_USER="${R_LIBS_USER:-$PROJECT_ROOT/.Rlib}"

# 1) R sceptre (canonical source)
"$RSCRIPT" "$ROOT/r_sceptre/run_sceptre_on_target_r.R" "$ROOT/r_sceptre"

# 2) Mirror R outputs into canonical sceptre/*.csv used by downstream scripts
ROOT_FOR_PY="$ROOT" python3 - <<'PY'
import pandas as pd
import os
from pathlib import Path
root=Path(os.environ["ROOT_FOR_PY"])
r = root/'r_sceptre'/'outputs'/'discovery_r_on_target.csv'
out = root/'sceptre'
out.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(r)
df.to_csv(out/'on_target_results.csv', index=False)
df.to_csv(out/'discovery_merged.csv', index=False)
pd.DataFrame(columns=df.columns).to_csv(out/'trans_panel_results.csv', index=False)
pd.DataFrame({'source':['R_sceptre'],'rows':[len(df)]}).to_csv(out/'result_source.csv', index=False)
print('Mirrored R outputs:', len(df), 'rows')
PY

# 3) Downstream analyses
python3 "$ROOT/run_on_target_expression_aware.py"
python3 "$ROOT/run_scanpy_wilcoxon_baseline.py"
python3 "$ROOT/run_scanpy_umap_leiden.py"
python3 "$ROOT/run_qc_marker_panels.py"
python3 "$ROOT/run_perturbation_burden.py"
python3 "$ROOT/run_pathway_hallmark.py"

echo "Completed R-canonical pipeline: $ROOT"
