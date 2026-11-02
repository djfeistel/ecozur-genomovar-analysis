#!/bin/bash

RESULTS="../results"

./plot_checkm_stats.py \
    -i ${RESULTS}/bin_stats_ext.tsv \
    -x 16 \
    -y 18 \
    -o "${RESULTS}/scaffolds_vs_genome_size.png" \
    --exclude ${RESULTS}/CheckM_failed_Completeness_ge90_Contamination_le5.txt \
    --highlight-excluded \
    --label-failed

