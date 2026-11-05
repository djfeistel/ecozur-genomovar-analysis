# EcoZUR genomovar assessment

## Files

- [EcoZUR.ani](inputs/EcoZUR.ani): pairwise fastANI comparisons. Columns are
  query genome, reference genome, ANI (%), matched fragments, and total query
  fragments.
- [EcoZUR_genomovar_metadata.tsv](inputs/EcoZUR_genomovar_metadata.tsv):
  genome identifiers, genomovar assignments based on 99.5% ANI and complete linkage, and phylogroups.
- [scripts/](scripts/): the analysis script, Bash launcher, and Python dependencies.
- [results/](results/): ANI violin plots, ANI and RBM heatmaps, numerical
  tables, and analysis reports. Figures are available as PNG, PDF, and SVG.

Raw reciprocal best match (RBM) files were generate using the [enveomics aai.rb](https://github.com/lmrodriguezr/enveomics/blob/main/Scripts/aai.rb) 
script. The files are not included here because the collection contains approximately **32,000 files**. The RBM heatmaps and summary tables
are included.

## How to run

With Python 3.12, run from this directory:

```bash
python3 -m pip install -r scripts/requirements.txt
bash scripts/run_genomovar_assessment.sh
```

The launcher runs the ANI analyses with seed 13, using paths relative to the
script. It overwrites ANI outputs in [results/](results/) but does not run RBM
analysis or modify existing RBM outputs.

To calculate RBM ratios, supply the external RBM directory:

```bash
python3 scripts/plot_within_genomovar_ani.py --rbm --rbm-dir /path/to/RBM
python3 scripts/plot_within_genomovar_ani.py --rbm-min-identity 99.8 --rbm-dir /path/to/RBM
```

The default RBM cutoff is 100%; `--rbm-min-identity` sets an inclusive cutoff.
Output names distinguish the two thresholds.

## Methods

**Genome selection:** Genomes without metadata are excluded from grouped
analyses. Genomovars with at least three members are shown, ordered by
phylogroup and then genomovar.

**ANI comparisons:** Self comparisons are removed. Each genome pair is
represented once, retaining the direction with the larger fifth-column
query-fragment count. Ties are selected randomly using seed 13.

**Violin plot:** Points show within-genomovar pairwise ANI, median bars show
the central value, and violins show the distribution. Labels give the number
of genomes (`n`).

**ANI heatmap:** Off-diagonal cells show mean ANI between members of two
genomovars. Diagonal cells show mean within-genomovar ANI.

**RBM heatmaps:** Each genome-pair ratio is the number of column-3 identities
at or above the cutoff divided by all RBM rows. Cells show the mean of these
ratios, giving each genome pair equal weight rather than pooling counts.
Diagonal cells summarize within-genomovar comparisons. Reciprocal files
count once, using the ANI-selected direction when both files are available.
Missing comparisons are excluded, not assigned zero.
