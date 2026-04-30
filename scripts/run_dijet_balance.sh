#!/usr/bin/env bash
# Run plots/softdrop/dijet_balance_check.py on every per-sample output
# under data-jets/. Outputs land in plots/softdrop/output/<sample>/.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="/Users/siddharthsingh/miniforge3/envs/pyRoot-env/bin/python"
SCRIPT="$REPO_ROOT/plots/softdrop/dijet_balance_check.py"

for d in "$REPO_ROOT/data-jets"/*/; do
    sample="$(basename "$d")"
    root_file=$(ls "$d"/alljets_*.root "$d"/dijets_*.root 2>/dev/null | head -1 || true)
    if [[ -z "$root_file" ]]; then
        echo "[skip] $sample: no jet ROOT"
        continue
    fi
    echo "========================================================================"
    echo "  $sample"
    echo "========================================================================"
    "$PY" "$SCRIPT" "$root_file"
    echo
done
