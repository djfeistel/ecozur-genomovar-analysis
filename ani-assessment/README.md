# ANI assessment

Scripts for checking how fastANI values behave across EcoZUR genome pairs and how well ANI thresholds recover MLST sequence types (STs).

```
ani-assessment/
├── input/     fastANI output and ST assignments
├── scripts/   analysis scripts and run wrappers
└── results/   figures, tables and reports written by the scripts
```

## Input

| File | Description |
|------|-------------|
| `input/EcoZUR.ani` | All-vs-all fastANI output for 248 genomes (tab-separated, no header): `query`, `reference`, `ANI`, `matched fragments`, `total query fragments`. Contains both directions (A vs B and B vs A) for every pair. |
| `input/EcoZUR_STs.tsv` | Two columns, `genome` and `ST`, with one row per genome. `XXX` marks an unknown/novel ST. |

## Pair filtering (shared method)

fastANI is not symmetric: A vs B and B vs A can give slightly different ANI values. Before any downstream analysis, each genome pair is reduced to a single ANI value:

1. Self comparisons (A vs A) are removed.
2. For each pair, keep the row with the larger total-query-fragment count (5th column). This means the larger genome is the query.
3. If both directions have the same fragment count, one is chosen at random using a fixed seed (`--seed`, default 13), so results are reproducible.
4. The two directions are **not** averaged; the kept row's ANI is used as is.

This is the same filtering used by `03.ANI/plot_within_genomovar_ani.py`.

## Scripts

### `ani_reciprocal_difference.py`

Measures how much fastANI differs between the two directions of each pair.

- Calculates the absolute difference between A vs B and B vs A for every pair and reports the mean, SD, median, min and max.
- `-p`: violin plot of the differences, both overall and split by the pair's mean ANI (bins set with `--ani-bins`).
- `-u`: scatter of the discarded direction against the kept direction (the larger-query rule above), with a trend line.
- `-o`: per-pair table of both ANI values and their difference.

Run with `bash scripts/run_ani_reciprocal_difference.sh`. Outputs are named `results/reciprocal_ani_differences.*`.

### `st_ani_f1_assessment.py`

Tests how well a single ANI cutoff separates genomes of the same ST from genomes of different STs. It uses only the filtered pairs described above.

- **Scoring (one ST at a time):** for a given ST, pairs where both genomes are that ST are "same-ST" pairs. Pairs where only one genome is that ST are "non-ST" pairs. Each pair is classified as TP, FP, TN or FN depending on whether its ANI is at or above the cutoff. Precision, recall, F1 and accuracy are then calculated from those counts (see the script docstring or report for the exact definitions). Genomes with ST `XXX` are only used as non-ST partners. The one genome without an ST (B37-6) is excluded.
- `--f1-plot`: makes a two-panel figure of F1 against the ANI cutoff, calculated every 0.1% ANI. Panel A covers 95–100% and panel B 98.8–100%. The 100% point is calculated at 99.9999%.
- Choosing which STs to plot:
  - `--top-sts N`: the N STs with the most genomes (default 10). If the cutoff falls in a tie, the lowest-numbered ST is kept and a warning is written to the report.
  - `--st-list FILE`: a single-column list of STs to plot.
- `--plot-pooled`: adds a line that combines all selected STs. It isn't drawn by default but is always included in the CSVs.
- `--vline [ANI]`: adds a dashed vertical reference line, at 99.5% if no value is given.
- `--csv-thresholds`: the cutoffs reported in the summary CSV (default 99, 99.2, 99.5, 99.8).

Run with `bash scripts/run_st_ani_f1_assessment.sh` (from inside `scripts/`), or call the script directly. Use `--help` for all options.

## Results

All outputs are written to `results/`:

| Output | From | Contents |
|--------|------|----------|
| `reciprocal_ani_differences.png` | `ani_reciprocal_difference.py` | Violin plot of direction differences |
| `reciprocal_ani_differences.tsv` | `ani_reciprocal_difference.py` | Per-pair ANI values in both directions |
| `reciprocal_ani_differences_summary.txt` | `ani_reciprocal_difference.py` | Summary statistics |
| `st_ani_f1_f1.{png,pdf,svg}` | `st_ani_f1_assessment.py` | Two-panel F1 vs ANI-cutoff figure |
| `st_ani_f1_metrics.csv` | `st_ani_f1_assessment.py` | One row per ST (plus pooled): genome count, number of pairs, and precision/recall/F1/accuracy at each CSV cutoff |
| `st_ani_f1_metrics_by_threshold.csv` | `st_ani_f1_assessment.py` | TP/FP/TN/FN and all metrics at every plotted cutoff |
| `st_ani_f1_report.txt` | `st_ani_f1_assessment.py` | Inputs, filtering counts, STs used, warnings and software versions |

## Requirements

Python 3 with numpy, pandas, matplotlib and seaborn.
