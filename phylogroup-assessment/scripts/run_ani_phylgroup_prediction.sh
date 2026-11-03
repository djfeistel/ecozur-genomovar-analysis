#!/usr/bin/env bash
set -euo pipefail

if (( $# != 0 )); then
    printf 'Usage: %s (no arguments; fixed seed 13)\n' "$0" >&2
    exit 2
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results"
PYTHON="${PYTHON:-python3}"

exec "$PYTHON" -B "${SCRIPT_DIR}/plot_clermontyping_ani.py" \
    --seed 13 \
    --output "${RESULTS_DIR}/ClermonTyping_ANI_clustermap.png" \
    --predictions-output "${RESULTS_DIR}/ClermonTyping_phylogroup_predictions.tsv"
