#!/usr/bin/env python3
"""Assess how well fastANI thresholds recover MLST sequence types (STs).

Run from any directory:
    python3 st_ani_f1_assessment.py --f1-plot                   # top 10 STs
    python3 st_ani_f1_assessment.py --f1-plot --top-sts 5
    python3 st_ani_f1_assessment.py --f1-plot --st-list my_sts.txt
    python3 st_ani_f1_assessment.py --f1-plot --vline           # line at 99.5%
    python3 st_ani_f1_assessment.py --f1-plot --vline 99.2
Requires numpy, pandas, matplotlib, and seaborn.
Defaults resolve relative to this script, not the current working directory.

Pair filtering (identical to 03.ANI/plot_within_genomovar_ani.py):
  Self comparisons are removed. For each unordered genome pair, the entire
  fastANI row with the largest fifth-column query-fragment count is retained
  (i.e. the larger genome is the query); ties are sampled uniformly using
  --seed. Reciprocal ANI values are NOT averaged. All metrics are computed
  only on these retained pairs.

ST assignment:
  STs come from the "ST" column of --metadata. Genomes in the ANI input that
  lack metadata are excluded (with a warning). Genomes whose ST is "XXX"
  (unknown/novel) are never evaluated as their own ST but are kept as
  non-ST partners (negatives) for the other STs.

Per-ST confusion counts (one-vs-rest), for ST X at ANI threshold t:
  Only retained pairs that contain at least one ST-X genome are used.
    same-ST pair : both genomes are ST X
    non-ST pair  : exactly one genome is ST X (the other is any other ST,
                   including "XXX")
    TP = same-ST pair and ANI >= t
    FP = same-ST pair and ANI <  t
    TN = non-ST pair  and ANI <  t
    FN = non-ST pair  and ANI >= t
  Precision = TP / (TP + FP)
  Recall    = TP / (TP + FN)
  F1        = 2 * Precision * Recall / (Precision + Recall)
            (computed as 2TP / (2TP + FP + FN), which is identical when
             both terms are defined and equals 0 when TP = 0)
  Accuracy  = (TP + TN) / (TP + TN + FP + FN)
  An undefined ratio (zero denominator) is reported as NA.

  The "All selected STs (pooled)" row sums TP, FP, TN and FN across the
  selected STs before computing the metrics. A non-ST pair linking two
  selected STs is counted once from each ST's side. It is always written to
  the CSVs but only drawn on the plot with --plot-pooled.

ST selection (only STs with >= 2 genomes have same-ST pairs and are eligible):
  --top-sts N    the N STs with the most genomes (ties broken naturally by ST)
  --st-list FILE single-column file of STs ("#" comments, blank lines and an
                 optional "ST" header are ignored)

--f1-plot outputs (prefix set with --output-prefix):
  <prefix>_f1.{png,pdf,svg}  two panels, one row: F1 vs ANI threshold over
                             95-100% (A) and 98.8-100% (B), step 0.1%. The
                             final threshold is 99.9999% (shown as 100 on
                             the axis) because no pair can exceed 100%.
                             --vline [ANI] adds a dashed black reference
                             line (99.5% if no value is given).
  <prefix>_metrics.csv       one row per ST (+ pooled): total genomes,
                             datapoints (TP+TN+FP+FN), same-ST and non-ST
                             pair counts, and precision/recall/F1/accuracy at
                             each --csv-thresholds value (default
                             99, 99.2, 99.5, 99.8).
  <prefix>_metrics_by_threshold.csv
                             long-format TP/FP/TN/FN and metrics for every
                             threshold in the plotting grid.
  <prefix>_report.txt        inputs, filtering counts, warnings, outputs.
"""

import argparse
from pathlib import Path
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MultipleLocator
import numpy as np
import pandas as pd
import seaborn as sns


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[2]
ANI_COLUMNS = ["query", "reference", "ani", "matched_fragments", "query_fragments"]
PAIR_COLUMNS = ["genome_a", "genome_b"]
UNKNOWN_ST = "XXX"
POOLED_LABEL = "All selected STs (pooled)"
MAX_THRESHOLD = 99.9999
PANELS = [(95.0, "A"), (98.8, "B")]
GRID_MIN = min(xmin for xmin, _ in PANELS)
GRID_STEP = 0.1
DEFAULT_CSV_THRESHOLDS = "99,99.2,99.5,99.8"
DEFAULT_VLINE = 99.5
METRICS = ["precision", "recall", "f1", "accuracy"]


def natural_key(label):
    return tuple(
        (1, int(part)) if part.isdigit() else (0, part.casefold())
        for part in re.split(r"(\d+)", str(label))
        if part
    )


# --- Filtering: copied unchanged from 03.ANI/plot_within_genomovar_ani.py ---

def read_ani(path):
    ani = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, keep_default_na=False)
    if ani.shape[1] != len(ANI_COLUMNS):
        raise ValueError(f"{path}: expected exactly five fastANI columns.")
    ani.columns = ANI_COLUMNS
    if ani.empty or ani[["query", "reference"]].eq("").any().any():
        raise ValueError(f"{path}: ANI input is empty or contains blank genome IDs.")
    for column in ANI_COLUMNS[2:]:
        ani[column] = pd.to_numeric(ani[column], errors="raise")
    if not np.isfinite(ani[ANI_COLUMNS[2:]].to_numpy()).all():
        raise ValueError(f"{path}: ANI and fragment counts must be finite.")
    if not ani["ani"].between(0, 100).all():
        raise ValueError(f"{path}: ANI must be a percentage between 0 and 100.")
    for column in ["matched_fragments", "query_fragments"]:
        values = ani[column]
        if (values < 0).any() or (values % 1 != 0).any():
            raise ValueError(f"{path}: {column} must contain nonnegative integers.")
        ani[column] = values.astype("int64")
    if (ani["query_fragments"] == 0).any():
        raise ValueError(f"{path}: query-fragment counts must be positive.")
    if (ani["matched_fragments"] > ani["query_fragments"]).any():
        raise ValueError(f"{path}: matched fragments exceed query fragments.")
    if ani.duplicated(["query", "reference"]).any():
        raise ValueError(f"{path}: duplicate directed genome comparisons.")
    return ani


def select_pairs(ani, seed):
    pairs = ani.loc[ani["query"] != ani["reference"]].copy()
    if pairs.empty:
        raise ValueError("No non-self ANI comparisons remain.")
    pairs[PAIR_COLUMNS] = np.sort(pairs[["query", "reference"]].to_numpy(), axis=1)
    maximum = pairs.groupby(PAIR_COLUMNS)["query_fragments"].transform("max")
    candidates = pairs.loc[pairs["query_fragments"] == maximum]
    grouped = candidates.groupby(PAIR_COLUMNS, sort=True)
    tied_pairs = int((grouped.size() > 1).sum())
    selected = grouped.sample(n=1, random_state=np.random.default_rng(seed))
    selected = selected.sort_values(PAIR_COLUMNS).reset_index(drop=True)
    return selected, tied_pairs


# --- ST annotation and selection ---

def read_st_metadata(path):
    metadata = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = {"genome", "ST"} - set(metadata.columns)
    if missing:
        raise ValueError(f"{path}: missing metadata columns: {', '.join(sorted(missing))}")
    metadata = metadata[["genome", "ST"]].apply(lambda column: column.str.strip())
    if metadata.empty or metadata.eq("").any().any():
        raise ValueError(f"{path}: genome or ST fields are empty.")
    if metadata["genome"].duplicated().any():
        raise ValueError(f"{path}: duplicate genome IDs in metadata.")
    return metadata


def annotate_pairs(selected, metadata, input_genomes):
    warnings = []
    st_of = metadata.set_index("genome")["ST"]
    no_metadata = sorted(input_genomes - set(st_of.index), key=natural_key)
    if no_metadata:
        warnings.append("Excluded (no metadata): " + ", ".join(no_metadata))
    pairs = selected.copy()
    pairs["st_a"] = pairs["genome_a"].map(st_of)
    pairs["st_b"] = pairs["genome_b"].map(st_of)
    excluded_pairs = int(pairs[["st_a", "st_b"]].isna().any(axis=1).sum())
    pairs = pairs.dropna(subset=["st_a", "st_b"]).reset_index(drop=True)
    genomes = metadata.loc[metadata["genome"].isin(input_genomes)]
    return pairs, genomes, excluded_pairs, warnings


def st_genome_counts(genomes):
    counts = genomes.loc[genomes["ST"] != UNKNOWN_ST, "ST"].value_counts()
    order = sorted(counts.index, key=lambda st: (-counts[st], natural_key(st)))
    return counts.reindex(order)


def read_st_list(path):
    sts = []
    for line in path.read_text().splitlines():
        value = line.split("#", 1)[0].strip()
        if value:
            sts.append(value.split()[0])
    if sts and sts[0].casefold() == "st":
        sts = sts[1:]
    sts = list(dict.fromkeys(sts))
    if not sts:
        raise ValueError(f"{path}: no STs listed.")
    return sts


def choose_sts(counts, top_sts, st_list):
    warnings = []
    eligible = counts[counts >= 2]
    if st_list is not None:
        requested = read_st_list(st_list)
        chosen = []
        for st in requested:
            if st == UNKNOWN_ST:
                warnings.append(f"ST {st} is unknown/novel and is never evaluated; skipped.")
            elif st not in counts.index:
                warnings.append(f"ST {st} is not present among the ANI genomes; skipped.")
            elif st not in eligible.index:
                warnings.append(f"ST {st} has 1 genome (no same-ST pairs); skipped.")
            else:
                chosen.append(st)
        if not chosen:
            raise ValueError(f"{st_list}: none of the listed STs can be evaluated.")
        return chosen, warnings
    if top_sts > len(eligible):
        warnings.append(
            f"Requested top {top_sts} STs but only {len(eligible)} have >= 2 genomes."
        )
    chosen = eligible.index[:top_sts].tolist()
    if len(eligible) > top_sts and eligible.iloc[top_sts] == eligible.iloc[top_sts - 1]:
        tied = eligible[eligible == eligible.iloc[top_sts - 1]].index.tolist()
        warnings.append(
            f"Top-{top_sts} cutoff falls within a tie at {eligible.iloc[top_sts - 1]} genomes "
            f"(STs {', '.join(tied)}); kept the first by natural ST order."
        )
    return chosen, warnings


# --- Metrics ---

def threshold_grid():
    steps = np.arange(round((100 - GRID_MIN) / GRID_STEP))
    return np.append(np.round(GRID_MIN + steps * GRID_STEP, 6), MAX_THRESHOLD)


def count_at_least(sorted_values, thresholds):
    return len(sorted_values) - np.searchsorted(sorted_values, thresholds, side="left")


def confusion_counts(pairs, st, thresholds):
    in_a = pairs["st_a"].eq(st)
    in_b = pairs["st_b"].eq(st)
    same = np.sort(pairs.loc[in_a & in_b, "ani"].to_numpy())
    cross = np.sort(pairs.loc[in_a ^ in_b, "ani"].to_numpy())
    tp = count_at_least(same, thresholds)
    fn = count_at_least(cross, thresholds)
    return pd.DataFrame({
        "threshold": thresholds,
        "TP": tp,
        "FP": len(same) - tp,
        "TN": len(cross) - fn,
        "FN": fn,
    })


def safe_divide(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denominator > 0, numerator / denominator, np.nan)


def add_metrics(table):
    tp, fp, tn, fn = (table[column] for column in ["TP", "FP", "TN", "FN"])
    table["precision"] = safe_divide(tp, tp + fp)
    table["recall"] = safe_divide(tp, tp + fn)
    table["f1"] = safe_divide(2 * tp, 2 * tp + fp + fn)
    table["accuracy"] = safe_divide(tp + tn, tp + tn + fp + fn)
    return table


def compute_metrics(pairs, sts, counts, thresholds):
    tables = []
    for st in sts:
        table = confusion_counts(pairs, st, thresholds)
        table.insert(0, "ST", st)
        table.insert(1, "total_genomes", int(counts[st]))
        tables.append(table)
    per_st = pd.concat(tables, ignore_index=True)
    pooled = per_st.groupby("threshold", sort=True)[["TP", "FP", "TN", "FN"]].sum()
    pooled = pooled.reset_index()
    pooled.insert(0, "ST", POOLED_LABEL)
    pooled.insert(1, "total_genomes", int(counts[sts].sum()))
    long = pd.concat([per_st, pooled], ignore_index=True)
    long.insert(2, "datapoints", long[["TP", "FP", "TN", "FN"]].sum(axis=1))
    long.insert(3, "same_st_pairs", long["TP"] + long["FP"])
    long.insert(4, "non_st_pairs", long["TN"] + long["FN"])
    return add_metrics(long)


def threshold_label(value):
    return np.format_float_positional(value, trim="-")


def summary_table(long, csv_thresholds):
    keys = ["ST", "total_genomes", "datapoints", "same_st_pairs", "non_st_pairs"]
    summary = long.drop_duplicates("ST")[keys].reset_index(drop=True)
    for threshold in csv_thresholds:
        rows = long.loc[long["threshold"] == threshold].set_index("ST")[METRICS]
        rows.columns = [f"{metric}_{threshold_label(threshold)}" for metric in METRICS]
        summary = summary.merge(rows, left_on="ST", right_index=True, how="left")
    return summary


# --- Plotting (ggplot2 theme_bw style) ---

def set_theme_bw():
    sns.set_theme(style="white", context="paper", font="DejaVu Sans")
    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.8,
        "axes.labelsize": 11,
        "axes.labelcolor": "black",
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.color": "#4D4D4D",
        "ytick.color": "#4D4D4D",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.fontsize": 8.5,
        "legend.title_fontsize": 9.5,
    })


def style_panel(ax, xmin):
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#333333")
        spine.set_linewidth(0.8)
    span = 100 - xmin
    major, minor = (0.2, 0.1) if span <= 2 else (1, 0.5) if span <= 5 else (2, 1)
    ax.set_xlim(xmin - 0.02 * span, 100 + 0.02 * span)
    ax.set_ylim(-0.03, 1.03)
    ax.xaxis.set_major_locator(MultipleLocator(major))
    ax.xaxis.set_minor_locator(MultipleLocator(minor))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(MultipleLocator(0.125))
    ax.grid(which="major", color="#EBEBEB", linewidth=0.7)
    ax.grid(which="minor", color="#F5F5F5", linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    ax.set_axisbelow(True)
    ax.set_xlabel("ANI threshold (%)")


def st_palette(n):
    if n <= 8:
        okabe_ito = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                     "#0072B2", "#D55E00", "#CC79A7", "#999999"]
        return okabe_ito[:n]
    if n <= 10:
        return sns.color_palette("colorblind", n)
    return sns.color_palette("husl", n)


def make_f1_figure(long, sts, counts, plot_pooled=False, vline=None):
    set_theme_bw()
    colors = dict(zip(sts, st_palette(len(sts))))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
    for ax, (xmin, tag) in zip(axes, PANELS):
        for st in sts:
            data = long.loc[long["ST"] == st]
            ax.plot(data["threshold"], data["f1"], color=colors[st], linewidth=1.2,
                    label=f"ST{st} (n = {counts[st]})", solid_capstyle="round", zorder=3)
        if plot_pooled:
            pooled = long.loc[long["ST"] == POOLED_LABEL]
            ax.plot(pooled["threshold"], pooled["f1"], color="black", linewidth=1.5,
                    linestyle=(0, (4, 2)), label=POOLED_LABEL, zorder=4)
        style_panel(ax, xmin)
        if vline is not None and xmin <= vline <= 100:
            ax.axvline(vline, color="black", linewidth=0.9, linestyle=(0, (5, 3)), zorder=2)
        ax.text(-0.02, 1.04, tag, transform=ax.transAxes, fontsize=14,
                fontweight="bold", ha="right", va="bottom")
        ax.set_title(f"{xmin:g}–100% ANI", fontsize=10.5, loc="left", pad=6)
    axes[0].set_ylabel("F1 score")
    handles, labels = axes[0].get_legend_handles_labels()
    if vline is not None:
        handles.append(Line2D([], [], color="black", linewidth=0.9, linestyle=(0, (5, 3))))
        labels.append(f"{vline:g}% ANI")
    fig.legend(handles, labels, title="Sequence type", loc="center left",
               bbox_to_anchor=(0.825, 0.5), frameon=False, handlelength=2.4,
               borderaxespad=0)
    fig.subplots_adjust(left=0.065, right=0.81, bottom=0.2, top=0.88, wspace=0.08)
    fig.text(0.065, 0.035,
             f"F1 evaluated every {GRID_STEP:g}% ANI; the 100% point is evaluated at "
             f"{MAX_THRESHOLD}% ANI. One retained fastANI value per genome pair "
             "(larger query-fragment count).",
             fontsize=7.5, color="#4D4D4D")
    return fig


# --- CLI ---

def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def nonnegative_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be a nonnegative integer")
    return number


def threshold_list(value):
    try:
        thresholds = [float(part) for part in value.split(",") if part.strip()]
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be comma-separated numbers") from error
    if not thresholds or not all(0 < t <= MAX_THRESHOLD for t in thresholds):
        raise argparse.ArgumentTypeError(f"thresholds must be in (0, {MAX_THRESHOLD}]")
    return thresholds


def ani_value(value):
    number = float(value)
    if not 0 < number <= 100:
        raise argparse.ArgumentTypeError("ANI value must be in (0, 100]")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ani", type=Path, default=SCRIPT_DIR.parent / "input" / "EcoZUR.ani",
                        help="Five-column fastANI input.")
    parser.add_argument("--metadata", type=Path,
                        default=SCRIPT_DIR.parent / "input" / "EcoZUR_STs.tsv",
                        help="TSV with 'genome' and 'ST' columns.")
    parser.add_argument("--output-prefix", type=Path,
                        default=SCRIPT_DIR.parent / "results" / "st_ani_f1",
                        help="Output path prefix (no extension; existing outputs overwritten).")
    parser.add_argument("--seed", type=nonnegative_int, default=13,
                        help="Random seed for tied query-fragment counts (default: 13).")
    parser.add_argument("--f1-plot", action="store_true",
                        help="Make the two-panel F1 plot and the per-ST metrics CSV files.")
    parser.add_argument("--plot-pooled", action="store_true",
                        help="Also draw the 'All selected STs (pooled)' F1 line "
                             "(default: not plotted; always written to the CSVs).")
    parser.add_argument("--vline", type=ani_value, nargs="?", const=DEFAULT_VLINE,
                        metavar="ANI",
                        help=f"Draw a dashed black vertical line at this ANI value "
                             f"(no value given = {DEFAULT_VLINE:g}%%; flag omitted = no line).")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--top-sts", type=positive_int, default=10,
                           help="Plot the N STs with the most genomes (default: 10).")
    selection.add_argument("--st-list", type=Path,
                           help="Single-column file of STs to plot (overrides --top-sts).")
    parser.add_argument("--csv-thresholds", type=threshold_list,
                        default=threshold_list(DEFAULT_CSV_THRESHOLDS),
                        help=f"Comma-separated ANI thresholds for the summary CSV "
                             f"(default: {DEFAULT_CSV_THRESHOLDS}).")
    args = parser.parse_args(argv)

    ani = read_ani(args.ani)
    metadata = read_st_metadata(args.metadata)
    input_genomes = set(ani["query"]) | set(ani["reference"])
    selected, tied_pairs = select_pairs(ani, args.seed)
    pairs, genomes, excluded_pairs, warnings = annotate_pairs(selected, metadata, input_genomes)
    counts = st_genome_counts(genomes)
    report = [
        f"ANI input: {args.ani.resolve()}",
        f"Metadata: {args.metadata.resolve()}",
        f"Seed: {args.seed}",
        f"Input ANI rows: {len(ani)}",
        f"Unique input genomes: {len(input_genomes)}",
        f"Self comparisons removed: {int(ani['query'].eq(ani['reference']).sum())}",
        f"Retained unordered pairs (larger query-fragment count): {len(selected)}",
        f"Equal-maximum fragment-count pairs randomized: {tied_pairs}",
        f"Pairs excluded (genome without metadata): {excluded_pairs}",
        f"Pairs available for ST analysis: {len(pairs)}",
        f"Genomes with ST '{UNKNOWN_ST}' (negatives only): "
        f"{int(genomes['ST'].eq(UNKNOWN_ST).sum())}",
        f"Known STs among ANI genomes: {len(counts)} "
        f"({int((counts >= 2).sum())} with >= 2 genomes)",
    ]
    if not args.f1_plot:
        for warning in warnings:
            print(f"WARNING: {warning}", file=sys.stderr)
        print("\n".join(report))
        print("\nNo plot requested; add --f1-plot to compute metrics and plot F1 curves.")
        return 0

    sts, selection_warnings = choose_sts(counts, args.top_sts, args.st_list)
    warnings.extend(selection_warnings)
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    thresholds = threshold_grid()
    long = compute_metrics(pairs, sts, counts, thresholds)
    csv_thresholds = args.csv_thresholds
    summary = summary_table(
        compute_metrics(pairs, sts, counts, np.array(csv_thresholds)), csv_thresholds
    )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    summary_path = Path(f"{args.output_prefix}_metrics.csv")
    summary.to_csv(summary_path, index=False, na_rep="NA", float_format="%.6f")
    outputs.append(summary_path)
    long_path = Path(f"{args.output_prefix}_metrics_by_threshold.csv")
    long.to_csv(long_path, index=False, na_rep="NA", float_format="%.6f")
    outputs.append(long_path)
    fig = make_f1_figure(long, sts, counts, plot_pooled=args.plot_pooled, vline=args.vline)
    for extension in ["png", "pdf", "svg"]:
        path = Path(f"{args.output_prefix}_f1.{extension}")
        fig.savefig(path, dpi=450, facecolor="white")
        outputs.append(path)
    plt.close(fig)

    selection_text = (f"ST list file {args.st_list.resolve()}" if args.st_list
                      else f"top {args.top_sts} STs by genome count")
    report += [
        f"ST selection: {selection_text}",
        "STs evaluated: " + ", ".join(f"ST{st} (n={counts[st]})" for st in sts),
        f"Plot grid: {GRID_MIN:g} to {thresholds[-2]:g}% ANI every {GRID_STEP:g}%, "
        f"plus {MAX_THRESHOLD}%",
        "CSV thresholds: " + ", ".join(threshold_label(t) for t in csv_thresholds),
        "",
        "Per ST X (one-vs-rest): only retained pairs containing an ST-X genome are used.",
        "  TP = both ST X, ANI >= t;  FP = both ST X, ANI < t",
        "  TN = one ST X, ANI < t;    FN = one ST X, ANI >= t",
        "Precision = TP/(TP+FP); Recall = TP/(TP+FN); F1 = 2PR/(P+R) = 2TP/(2TP+FP+FN);",
        "Accuracy = (TP+TN)/(TP+TN+FP+FN). Zero denominators give NA.",
        f"Genomes with ST '{UNKNOWN_ST}' are negatives only, never evaluated as an ST.",
        "Pooled: TP/FP/TN/FN summed across the selected STs; non-ST pairs between two",
        "selected STs are counted once from each side.",
        "",
        *[f"WARNING: {warning}" for warning in warnings],
        "",
        f"Versions: numpy={np.__version__}, pandas={pd.__version__}, "
        f"matplotlib={matplotlib.__version__}, seaborn={sns.__version__}",
    ]
    report_path = Path(f"{args.output_prefix}_report.txt")
    outputs.append(report_path)
    report += ["", "Outputs:", *[str(path.resolve()) for path in outputs]]
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
