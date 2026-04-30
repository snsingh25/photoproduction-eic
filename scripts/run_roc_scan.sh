#!/usr/bin/env bash
# Run plots/softdrop/roc_eta_scan.py on every per-sample output under
# data-jets/. The script auto-routes outputs (per-eta PDFs, grid,
# roc_summary.log) into plots/softdrop/output/<sample>/.
#
# Usage:  ./scripts/run_roc_scan.sh
#         ./scripts/run_roc_scan.sh --bins "-1,0,1,2,3"
#
# Prerequisites: pyRoot-env conda env (has numpy + uproot + matplotlib).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="/Users/siddharthsingh/miniforge3/envs/pyRoot-env/bin/python"
SCRIPT="$REPO_ROOT/plots/softdrop/roc_eta_scan.py"

EXTRA_ARGS=("$@")

for d in "$REPO_ROOT/data-jets"/*/; do
    sample="$(basename "$d")"
    root_file=$(ls "$d"/alljets_*.root "$d"/dijets_*.root 2>/dev/null | head -1 || true)
    if [[ -z "$root_file" ]]; then
        echo "[skip] $sample: no jet ROOT found"
        continue
    fi
    echo "========================================================================"
    echo "  $sample"
    echo "========================================================================"
    "$PY" "$SCRIPT" "$root_file" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
    echo
done

echo "Done. Outputs under plots/softdrop/output/<sample>/."
