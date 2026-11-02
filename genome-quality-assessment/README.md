# Genome quality assessment (CheckM)

Quality assessment of the EcoZUR *E. coli* genome assemblies using CheckM
completeness, contamination and assembly statistics.

## Contents

```
genome-quality-assessment/
├── results/
│   ├── bin_stats_ext.tsv                                      # CheckM statistics for 265 genomes
│   ├── CheckM_failed_Completeness_ge90_Contamination_le5.txt  # 10 genomes failing the quality thresholds
│   ├── contamination_vs_completeness.png
│   └── scaffolds_vs_genome_size.png
└── scripts/
    ├── 01_run_completeness_vs_contamination.sh
    ├── 02_run_scaffolds_vs_genome_size.sh
    └── plot_checkm_stats.py                                   # plotting script (run with -h for full help)
```

## Data

`results/bin_stats_ext.tsv` is the official CheckM output for **265 genomes** from
the EcoZUR study. This set excludes the 14 genomes that did not assemble and was
generated **before** dereplication of genomes from the same host at >99.9% ANI.

`results/CheckM_failed_Completeness_ge90_Contamination_le5.txt` lists the **10
genomes** that did not meet the quality thresholds of >90% completeness and <5%
contamination. It was generated with the `--write-failed` flag of
`plot_checkm_stats.py`:

```bash
./plot_checkm_stats.py -i ../results/bin_stats_ext.tsv \
    -x Completeness -y Contamination --xline 90 --yline 5 \
    --write-failed ../results/CheckM_failed_Completeness_ge90_Contamination_le5.txt
```

This file is then passed to the plotting scripts with
`--exclude ... --highlight-excluded --label-failed`, so these 10 genomes are
shown as labelled orange diamonds in every plot.

## Usage

Run the scripts in order of their `01_` / `02_` prefix, from inside the
`scripts/` directory (they use relative paths to `../results`):

```bash
cd scripts
./01_run_completeness_vs_contamination.sh   # -> results/contamination_vs_completeness.png
./02_run_scaffolds_vs_genome_size.sh        # -> results/scaffolds_vs_genome_size.png
```

Requires Python 3 with `pandas`, `numpy` and `matplotlib`.
