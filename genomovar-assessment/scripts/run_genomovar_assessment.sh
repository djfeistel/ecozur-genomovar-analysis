#!/usr/bin/env bash
set -euo pipefail

if (( $# != 0 )); then
    printf 'Usage: %s (no arguments; ANI analyses only, seed 13)\n' "$0" >&2
    exit 2
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"

exec "$PYTHON" -B "${SCRIPT_DIR}/plot_within_genomovar_ani.py" \
    --ani "${SCRIPT_DIR}/../inputs/EcoZUR.ani" \
    --metadata "${SCRIPT_DIR}/../inputs/EcoZUR_genomovar_metadata.tsv" \
    --output-prefix "${SCRIPT_DIR}/../results/within_genomovar_ANI" \
    --seed 13
