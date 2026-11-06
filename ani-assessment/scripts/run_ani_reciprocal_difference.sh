#!/bin/bash
# Reproduce the reciprocal fastANI difference results:
#   - summary statistics (mean, SD, median, min, max of |ANI(A vs B) - ANI(B vs A)|)
#   - violin plot (all pairs + mean-ANI bins 97, 98, 99, 99.5), without point overlay
# Usage: bash run_ani_reciprocal_difference.sh   (can be run from any directory)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

INPUT="../input/EcoZUR.ani"
RESULTS="../results"
ANI_BINS="97,98,99,99.5"

mkdir -p "${RESULTS}"

python3 ani_reciprocal_difference.py \
    -i "${INPUT}" \
    -p "${RESULTS}/reciprocal_ani_differences.png" \
    --ani-bins "${ANI_BINS}" \
    -o ${RESULTS}/reciprocal_ani_differences.tsv \
    | tee "${RESULTS}/reciprocal_ani_differences_summary.txt"
