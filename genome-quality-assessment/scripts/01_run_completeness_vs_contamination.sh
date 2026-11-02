#!/bin/bash

RESULTS="../results"

./plot_checkm_stats.py \
    -i ${RESULTS}/bin_stats_ext.tsv \
    -x 12 \
    -y 13 \
    -o "${RESULTS}/contamination_vs_completeness.png" \
    --xline 90 \
    --yline 5 \
    --xlim 88 100 \
    --ylim 0 80 \
    --pad 0.02 \
    --label-failed \
    --exclude ${RESULTS}/CheckM_failed_Completeness_ge90_Contamination_le5.txt \
    --highlight-excluded 
