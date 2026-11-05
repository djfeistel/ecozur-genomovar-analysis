#!/usr/bin/env python3
"""Plot within-genomovar fastANI distributions, ordered by phylogroup.

Run from any directory:
    bash /path/to/genomovar-assessment/scripts/run_genomovar_assessment.sh

Requires numpy, pandas, matplotlib, and seaborn >= 0.13.
Defaults resolve relative to this script, not the current working directory.
Bundled ANI and metadata are in ../inputs/; outputs are written to ../results/.
Use --help for input/output options.

Self comparisons are removed. For each unordered genome pair, retain the
entire fastANI row with the largest fifth-column query-fragment count, NOT
the largest ANI. Ties are sampled uniformly using --seed (default 13).
The selected ANI values are not averaged across reciprocal directions.

Only annotated genomes present in the ANI input count toward the minimum
of three members per genomovar. Unannotated genomes are reported and excluded
from grouped analysis, but remain in the selected-pairs export. Genomovars
must map to a single phylogroup. Both labels are naturally sorted (g9 < g10).

Outputs: all selected pairs, plotted within-genomovar pairs, descriptive
statistics, a text analysis report, and PNG (450 dpi), PDF, and SVG figures.
The seed also controls horizontal point jitter. Each point is one unordered
pair, not an independent biological replicate; violins are descriptive KDEs
clipped to the observed range, with equal maximum widths and median bars.

A separate genomovar matrix section uses the same selected genome pairs and
the exact violin-axis order. Off-diagonal cells are unweighted mean ANI across
members of two different genomovars; diagonal cells are the within-genomovar
means, not 100%. Additional outputs use the "_genomovar_matrix" suffix: a mean
ANI TSV, pair-count TSV, method report, and annotated PNG/PDF/SVG heatmaps.

Optional RBM analysis:
    python3 scripts/plot_within_genomovar_ani.py --rbm --rbm-dir /path/to/RBM
    python3 scripts/plot_within_genomovar_ani.py --rbm-min-identity 99.8 --rbm-dir /path/to/RBM

--rbm uses a 100% cutoff; supplying --rbm-min-identity also enables this section.
Raw RBM files are not bundled; supply their external location with --rbm-dir.
Each file's ratio is the number of column-3 values >= the cutoff divided by
its total number of nonblank RBM rows. Genomovar cells average genome-pair
ratios with equal weight, not pooled RBM counts. Reciprocal files count once:
use the retained ANI direction if both exist, otherwise the available file.
Missing comparisons are reported, never treated as zero. Only genomes in the
existing violin analysis are used; genome self comparisons are excluded.
RBM files are streamed with bounded memory. Threshold-specific output names
keep different cutoffs separate. Identity is used as reported in column 3;
100% reported identity does not establish full-length gene identity.
"""

import argparse
from itertools import combinations
from pathlib import Path
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import numpy as np
import pandas as pd
import seaborn as sns


SCRIPT_DIR = Path(__file__).resolve().parent
ANI_COLUMNS = ["query", "reference", "ani", "matched_fragments", "query_fragments"]
METADATA_COLUMNS = ["genome", "genomovar_complete", "phylogroup"]
PAIR_COLUMNS = ["genome_a", "genome_b"]
GROUP_COLUMNS = ["phylogroup", "genomovar_complete"]
MIN_MEMBERS = 3


def natural_key(label):
    return tuple(
        (1, int(part)) if part.isdigit() else (0, part.casefold())
        for part in re.split(r"(\d+)", str(label))
        if part
    )


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


def read_metadata(path):
    metadata = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = set(METADATA_COLUMNS) - set(metadata.columns)
    if missing:
        raise ValueError(f"{path}: missing metadata columns: {', '.join(sorted(missing))}")
    metadata = metadata[METADATA_COLUMNS].copy()
    for column in METADATA_COLUMNS:
        metadata[column] = metadata[column].str.strip()
    if metadata.empty or metadata.eq("").any().any():
        raise ValueError(f"{path}: required metadata fields are empty.")
    if metadata["genome"].duplicated().any():
        raise ValueError(f"{path}: duplicate genome IDs in metadata.")
    return metadata


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


def summarize_groups(selected, metadata, input_genomes):
    warnings = []
    missing = sorted(input_genomes - set(metadata["genome"]), key=natural_key)
    absent = sorted(set(metadata["genome"]) - input_genomes, key=natural_key)
    if missing:
        warnings.append(
            "Excluded from grouped analysis (no metadata): " + ", ".join(missing)
        )
    if absent:
        warnings.append(
            "Not counted toward membership (absent from ANI): " + ", ".join(absent)
        )
    present = metadata.loc[metadata["genome"].isin(input_genomes)].copy()
    conflicts = present.groupby("genomovar_complete")["phylogroup"].nunique()
    conflicts = conflicts[conflicts > 1]
    if not conflicts.empty:
        raise ValueError(
            "Genomovars assigned to multiple phylogroups: "
            + ", ".join(conflicts.index)
        )
    membership = present.groupby(GROUP_COLUMNS).size().rename("n_genomes").reset_index()
    membership = membership.loc[membership["n_genomes"] >= MIN_MEMBERS].copy()
    if membership.empty:
        raise ValueError("No annotated genomovars have at least three input genomes.")
    order = sorted(
        membership.index,
        key=lambda i: (
            natural_key(membership.at[i, "phylogroup"]),
            natural_key(membership.at[i, "genomovar_complete"]),
        ),
    )
    membership = membership.loc[order].reset_index(drop=True)
    annotation = present.set_index("genome")
    query_group = selected["query"].map(annotation["genomovar_complete"])
    reference_group = selected["reference"].map(annotation["genomovar_complete"])
    eligible = membership["genomovar_complete"].tolist()
    within = selected.loc[
        query_group.eq(reference_group) & query_group.isin(eligible)
    ].copy()
    within["genomovar_complete"] = query_group.loc[within.index]
    within["phylogroup"] = within["query"].map(annotation["phylogroup"])
    rank = {label: i for i, label in enumerate(eligible)}
    within = within.sort_values(
        "genomovar_complete", key=lambda labels: labels.map(rank), kind="stable"
    ).reset_index(drop=True)
    stats = within.groupby(GROUP_COLUMNS)["ani"].agg(
        n_pairs="size",
        mean_ani="mean",
        median_ani="median",
        sd_ani="std",
        min_ani="min",
        q25_ani=lambda values: values.quantile(0.25),
        q75_ani=lambda values: values.quantile(0.75),
        max_ani="max",
    )
    summary = membership.merge(stats, on=GROUP_COLUMNS, how="left", validate="one_to_one")
    if summary["n_pairs"].isna().any():
        empty = summary.loc[summary["n_pairs"].isna(), "genomovar_complete"]
        raise ValueError("No within-genomovar ANI comparisons for: " + ", ".join(empty))
    summary["n_pairs"] = summary["n_pairs"].astype(int)
    summary.insert(
        4, "expected_pairs", summary["n_genomes"] * (summary["n_genomes"] - 1) // 2
    )
    incomplete = summary.loc[summary["n_pairs"] != summary["expected_pairs"]]
    for row in incomplete.itertuples(index=False):
        warnings.append(
            f"{row.genomovar_complete}: only {row.n_pairs}/{row.expected_pairs} "
            "within-genomovar pairs are available; no missing ANI values imputed."
        )
    return within, summary, warnings


def make_figure(within, summary, seed):
    sns.set_theme(style="whitegrid", context="paper", font="DejaVu Sans")
    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.labelsize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 10,
    })
    order = summary["genomovar_complete"].tolist()
    fig, ax = plt.subplots(figsize=(max(9, 0.44 * len(order) + 1.5), 6.4))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.25, top=0.81)
    sns.violinplot(
        data=within,
        x="genomovar_complete",
        y="ani",
        order=order,
        color="#D9D9D9",
        linecolor="#888888",
        linewidth=0.8,
        inner=None,
        cut=0,
        density_norm="width",
        bw_method="scott",
        width=0.8,
        ax=ax,
    )
    rng = np.random.default_rng(seed)
    for i, group in enumerate(order):
        values = within.loc[within["genomovar_complete"] == group, "ani"].to_numpy()
        ax.scatter(
            i + rng.uniform(-0.19, 0.19, len(values)),
            values,
            s=13,
            color="#444444",
            alpha=0.75,
            edgecolors="none",
            zorder=3,
        )
        ax.hlines(np.median(values), i - 0.23, i + 0.23, color="#222222",
                  linewidth=1.3, zorder=4)
    lower = float(within["ani"].min())
    upper = float(within["ani"].max())
    padding = max((upper - lower) * 0.06, 0.01)
    ax.set_ylim(lower - padding, min(100, upper + padding) + 0.005)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([
        f"{row.genomovar_complete}\n(n={row.n_genomes})"
        for row in summary.itertuples(index=False)
    ])
    ax.set_xlabel("Genomovar", labelpad=12)
    ax.set_ylabel("Pairwise ANI (%)", labelpad=10)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    if upper - lower < 1:
        ax.yaxis.set_major_locator(MultipleLocator(0.1 if upper - lower > 0.2 else 0.02))
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, color="#E5E5E5", linewidth=0.7)
    ax.set_axisbelow(True)
    sns.despine(ax=ax)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_color("#888888")
    for phylogroup, block in summary.groupby("phylogroup", sort=False):
        start, stop = int(block.index.min()), int(block.index.max())
        if start:
            ax.axvline(start - 0.5, color="#B5B5B5", linewidth=0.8, zorder=1)
        strip = Rectangle(
            (start - 0.5, 1.025), stop - start + 1, 0.08,
            transform=ax.get_xaxis_transform(), clip_on=False,
            facecolor="#EFEFEF", edgecolor="white", linewidth=1.5,
        )
        ax.add_patch(strip)
        ax.text(
            (start + stop) / 2, 1.065, f"Phylogroup {phylogroup}",
            transform=ax.get_xaxis_transform(), ha="center", va="center",
            fontsize=11, fontweight="semibold",
        )
    fig.suptitle("Within-genomovar average nucleotide identity",
                 fontsize=15, fontweight="semibold", y=0.96)
    fig.text(0.07, 0.085,
             "Points: unique genome pairs; bars: medians; violins: observed-range densities.",
             fontsize=9, color="#444444")
    fig.text(0.07, 0.05,
             "Genomovars with at least 3 genomes; reciprocal direction selected by larger "
             f"query-fragment count (random ties; seed {seed}).",
             fontsize=9, color="#444444")
    return fig


def nonnegative_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be a nonnegative integer")
    return number


def build_genomovar_matrix(selected, metadata, summary):
    order = summary["genomovar_complete"].tolist()
    rank = {label: i for i, label in enumerate(order)}
    annotation = metadata.set_index("genome")["genomovar_complete"]
    query_group = selected["query"].map(annotation)
    reference_group = selected["reference"].map(annotation)
    keep = query_group.isin(order) & reference_group.isin(order)
    pairs = selected.loc[keep, ["ani"]].copy()
    query_rank = query_group.loc[keep].map(rank).to_numpy(dtype=int)
    reference_rank = reference_group.loc[keep].map(rank).to_numpy(dtype=int)
    pairs["row"] = np.minimum(query_rank, reference_rank)
    pairs["column"] = np.maximum(query_rank, reference_rank)
    grouped = pairs.groupby(["row", "column"])["ani"].agg(mean_ani="mean", n_pairs="size")

    labels = pd.Index(order, name="genomovar_complete")
    matrix = pd.DataFrame(np.nan, index=labels, columns=labels)
    counts = pd.DataFrame(0, index=labels, columns=labels, dtype=int)
    for (i, j), row in grouped.iterrows():
        matrix.iat[i, j] = matrix.iat[j, i] = row["mean_ani"]
        counts.iat[i, j] = counts.iat[j, i] = int(row["n_pairs"])
    if not np.allclose(np.diag(matrix), summary["mean_ani"], rtol=0, atol=1e-10):
        raise ValueError("Matrix diagonal does not match the within-genomovar ANI summary.")

    warnings = []
    sizes = summary["n_genomes"].to_numpy()
    for i, group_a in enumerate(order):
        for j in range(i, len(order)):
            expected = sizes[i] * (sizes[i] - 1) // 2 if i == j else sizes[i] * sizes[j]
            observed = counts.iat[i, j]
            if observed > expected:
                raise ValueError(f"{group_a}/{order[j]}: more genome pairs than expected.")
            if observed < expected:
                detail = "mean is NA" if observed == 0 else "mean uses available pairs only"
                warnings.append(
                    f"{group_a}/{order[j]}: {observed}/{expected} genome pairs; "
                    f"{detail}; no ANI values imputed."
                )
    return matrix, counts, warnings


def make_matrix_figure(
    matrix, summary, *,
    title="Mean ANI within and between genomovars",
    colorbar_label="Mean pairwise ANI (%)",
    vmin=None, vmax=100, footnotes=None,
):
    sns.set_theme(style="white", context="paper", font="DejaVu Sans")
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    size = max(10, 0.27 * len(summary) + 4.5)
    fig = plt.figure(figsize=(size, size))
    ax = fig.add_axes((0.20, 0.19, 0.70, 0.70))
    colorbar_ax = fig.add_axes((0.92, 0.19, 0.018, 0.70))
    labels = [
        f"{row.genomovar_complete} (n={row.n_genomes})"
        for row in summary.itertuples(index=False)
    ]
    minimum = np.floor(np.nanmin(matrix.to_numpy()) * 10) / 10 if vmin is None else vmin
    sns.heatmap(
        matrix,
        ax=ax,
        cbar_ax=colorbar_ax,
        cmap="viridis",
        vmin=minimum,
        vmax=vmax,
        square=True,
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 6},
        linewidths=0.25,
        linecolor="#FFFFFF",
        xticklabels=labels,
        yticklabels=labels,
        cbar_kws={"label": colorbar_label},
    )
    ax.set_facecolor("#EEEEEE")
    for i, j in np.argwhere(matrix.isna().to_numpy()):
        ax.text(j + 0.5, i + 0.5, "NA", ha="center", va="center",
                fontsize=6, color="#555555")
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_yticklabels(labels, rotation=0, fontsize=8)
    ax.tick_params(length=0)
    ax.set_xlabel("Genomovar (n = genomes)", fontsize=12, labelpad=12)
    ax.set_ylabel("Genomovar", fontsize=12)
    ax.yaxis.set_label_coords(-0.25, 0.5)
    colorbar_ax.tick_params(labelsize=9)
    colorbar_ax.yaxis.label.set_size(11)
    for i in range(len(summary)):
        ax.add_patch(Rectangle((i, i), 1, 1, fill=False, edgecolor="#333333", linewidth=0.7))
    for phylogroup, block in summary.reset_index(drop=True).groupby("phylogroup", sort=False):
        start, stop = int(block.index.min()), int(block.index.max()) + 1
        if start:
            ax.axvline(start, color="white", linewidth=1.6)
            ax.axhline(start, color="white", linewidth=1.6)
        ax.add_patch(Rectangle(
            (start, 1.015), stop - start, 0.035,
            transform=ax.get_xaxis_transform(), clip_on=False,
            facecolor="#EFEFEF", edgecolor="white", linewidth=1.5,
        ))
        ax.text(
            (start + stop) / 2, 1.0325, f"Phylogroup {phylogroup}",
            transform=ax.get_xaxis_transform(), ha="center", va="center", fontsize=10,
        )
        ax.add_patch(Rectangle(
            (-0.19, start), 0.035, stop - start,
            transform=ax.get_yaxis_transform(), clip_on=False,
            facecolor="#EFEFEF", edgecolor="white", linewidth=1.5,
        ))
        ax.text(
            -0.1725, (start + stop) / 2, f"Phylogroup {phylogroup}",
            transform=ax.get_yaxis_transform(), ha="center", va="center",
            rotation=90, fontsize=10,
        )
    fig.suptitle(title, fontsize=16, fontweight="semibold", y=0.97)
    if footnotes is None:
        footnotes = (
            "Cells: mean ANI across unique genome pairs; outlined diagonal: within-genomovar means.",
            "Same genomovars, order, and retained comparisons as the violin plot; "
            "cell labels rounded to 2 decimals.",
        )
    for y, text in zip((0.055, 0.035), footnotes):
        fig.text(0.20, y, text, fontsize=9, color="#444444")
    return fig


def write_genomovar_matrix(selected, metadata, summary, output_prefix, seed):
    matrix, counts, warnings = build_genomovar_matrix(selected, metadata, summary)
    for warning in warnings:
        print(f"WARNING: Matrix: {warning}", file=sys.stderr)
    prefix = Path(f"{output_prefix}_genomovar_matrix")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    matrix_path = Path(f"{prefix}.tsv")
    counts_path = Path(f"{prefix}_pair_counts.tsv")
    matrix.to_csv(matrix_path, sep="\t", na_rep="NA")
    counts.to_csv(counts_path, sep="\t")
    outputs = [matrix_path, counts_path]
    fig = make_matrix_figure(matrix, summary)
    for extension in ["png", "pdf", "svg"]:
        output = Path(f"{prefix}.{extension}")
        fig.savefig(output, dpi=450, facecolor="white")
        outputs.append(output)
    plt.close(fig)
    report_path = Path(f"{prefix}_report.txt")
    outputs.append(report_path)
    report_lines = [
        "Genomovar ANI matrix",
        f"Seed used for the shared retained genome pairs: {seed}",
        f"Genomovars: {len(summary)}",
        f"Genomes: {int(summary['n_genomes'].sum())}",
        f"Unique genome pairs used: {int(np.triu(counts.to_numpy()).sum())}",
        f"Within-genomovar pairs: {int(np.trace(counts.to_numpy()))}",
        f"Between-genomovar pairs: {int(np.triu(counts.to_numpy(), k=1).sum())}",
        "Axis order: " + ", ".join(summary["genomovar_complete"]),
        "Phylogroup order: " + ", ".join(summary["phylogroup"].drop_duplicates()),
        "",
        "Off-diagonal: unweighted arithmetic mean of retained ANI values for all "
        "genome pairs spanning the two genomovars; expected count = n_i * n_j.",
        "Diagonal: unweighted within-genomovar mean, matching the violin summary; "
        "expected count = n_i * (n_i - 1) / 2. Self comparisons are excluded.",
        "The same selected genome-pair rows feed both analyses; no new reciprocal "
        "selection or randomization occurs in this section.",
        "Both matrices are symmetric. Mirrored entries represent the same observations, "
        "not additional independent pairs.",
        "Only genomovars included in the violin plot are used, in the same order.",
        "The TSV retains unrounded means; only heatmap annotations use two decimals.",
        "Missing comparisons are not imputed. An empty cell has ANI=NA and pair count=0.",
        f"Input paths, exclusions, and software versions: {output_prefix}_report.txt",
        "",
        *[f"WARNING: {warning}" for warning in warnings],
        "Outputs:",
        *[str(path.resolve()) for path in outputs],
    ]
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("\n" + "\n".join(report_lines))


def identity_percentage(value):
    number = float(value)
    if not 0 <= number <= 100:
        raise argparse.ArgumentTypeError("identity cutoff must be finite and between 0 and 100")
    return number


def count_rbm_identities(path, threshold):
    qualifying = total = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split(None, 3)
            if not fields:
                continue
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least three RBM columns.")
            try:
                identity = float(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: nonnumeric column-3 identity {fields[2]!r}."
                ) from error
            if not 0 <= identity <= 100:
                raise ValueError(
                    f"{path}:{line_number}: identity must be finite and between 0 and 100."
                )
            total += 1
            qualifying += identity >= threshold
    if total == 0:
        raise ValueError(f"{path}: no RBM rows; a ratio with zero denominator is undefined.")
    return qualifying, total


def read_rbm_pairs(rbm_dir, selected, metadata, summary, input_genomes, threshold):
    if not rbm_dir.is_dir():
        raise NotADirectoryError(f"RBM directory does not exist: {rbm_dir}")
    order = summary["genomovar_complete"].tolist()
    rank = {label: i for i, label in enumerate(order)}
    included = metadata.loc[
        metadata["genome"].isin(input_genomes) & metadata["genomovar_complete"].isin(order)
    ]
    group_of = included.set_index("genome")["genomovar_complete"].to_dict()
    preferred = {
        tuple(sorted((row.query, row.reference))): f"{row.query}___{row.reference}.rbm"
        for row in selected.itertuples(index=False)
    }
    files = {}
    inventory = {"files_found": 0, "self_files_excluded": 0, "out_of_scope_files_excluded": 0}
    for path in sorted(rbm_dir.glob("*.rbm")):
        inventory["files_found"] += 1
        parts = path.stem.split("___")
        if len(parts) != 2 or not all(parts):
            raise ValueError(f"{path}: expected filename genomeA___genomeB.rbm.")
        a, b = parts
        if a == b:
            inventory["self_files_excluded"] += 1
        elif a not in group_of or b not in group_of:
            inventory["out_of_scope_files_excluded"] += 1
        else:
            files.setdefault(tuple(sorted((a, b))), []).append(path)
    if not files:
        raise ValueError(f"{rbm_dir}: no RBM files for the plotted genomes.")
    inventory["reciprocal_pairs"] = sum(len(paths) == 2 for paths in files.values())
    warnings = []
    records = []
    processed = 0
    print(f"RBM: streaming {len(files):,} files at identity >= {threshold:g}%.",
          file=sys.stderr, flush=True)
    for a, b in combinations(sorted(group_of), 2):
        group_a, group_b = sorted((group_of[a], group_of[b]), key=rank.__getitem__)
        record = {
            "genome_a": a, "genome_b": b,
            "genomovar_a": group_a, "genomovar_b": group_b,
            "identity_cutoff": threshold,
            "rbm_file": "", "unused_reciprocal_file": "",
            "qualifying_rows": None, "total_rows": None, "ratio": np.nan,
            "status": "missing",
        }
        paths = files.get((a, b), [])
        if not paths:
            warnings.append(f"Missing RBM comparison: {a}/{b} ({group_a}/{group_b}).")
        else:
            if len(paths) == 1:
                path = paths[0]
            else:
                matches = [path for path in paths if path.name == preferred.get((a, b))]
                if len(matches) != 1:
                    raise ValueError(
                        f"{a}/{b}: reciprocal RBM files have no unique retained ANI direction."
                    )
                path = matches[0]
                record["unused_reciprocal_file"] = next(
                    other.name for other in paths if other != path
                )
            qualifying, total = count_rbm_identities(path, threshold)
            record.update(
                rbm_file=path.name, qualifying_rows=qualifying, total_rows=total,
                ratio=qualifying / total, status="available",
            )
            processed += 1
            if processed % 1000 == 0 or processed == len(files):
                print(f"RBM: processed {processed:,}/{len(files):,} files.",
                      file=sys.stderr, flush=True)
        records.append(record)
    pairs = pd.DataFrame(records)
    for column in ["qualifying_rows", "total_rows"]:
        pairs[column] = pairs[column].astype("Int64")
    return pairs, inventory, warnings


def build_rbm_matrix(pairs, summary):
    order = summary["genomovar_complete"].tolist()
    grouped = pairs.groupby(["genomovar_a", "genomovar_b"])["ratio"].agg(
        mean_ratio="mean", n_pairs="count"
    )
    means = grouped["mean_ratio"].unstack()
    counts = grouped["n_pairs"].unstack()
    matrix = means.combine_first(means.T).reindex(index=order, columns=order)
    counts = counts.combine_first(counts.T).reindex(index=order, columns=order).fillna(0).astype(int)
    for table in (matrix, counts):
        table.index.name = table.columns.name = "genomovar_complete"
    sizes = summary["n_genomes"].to_numpy()
    expected = np.outer(sizes, sizes)
    np.fill_diagonal(expected, sizes * (sizes - 1) // 2)
    if (counts.to_numpy() > expected).any():
        raise ValueError("RBM genome-pair counts exceed the expected genomovar membership.")
    warnings = []
    for i, j in np.argwhere(np.triu(counts.to_numpy() < expected)):
        observed = counts.iat[i, j]
        detail = "mean is NA" if observed == 0 else "mean uses available pairs only"
        warnings.append(
            f"{order[i]}/{order[j]}: {observed}/{expected[i, j]} RBM genome pairs; "
            f"{detail}; no missing comparisons treated as zero."
        )
    return matrix, counts, warnings


def write_rbm_matrix(selected, metadata, summary, input_genomes,
                     rbm_dir, output_prefix, threshold, seed):
    pairs, inventory, warnings = read_rbm_pairs(
        rbm_dir, selected, metadata, summary, input_genomes, threshold
    )
    matrix, counts, cell_warnings = build_rbm_matrix(pairs, summary)
    warnings.extend(cell_warnings)
    for warning in warnings:
        print(f"WARNING: RBM: {warning}", file=sys.stderr)
    cutoff_label = np.format_float_positional(threshold, trim="-")
    prefix = Path(f"{output_prefix}_rbm_identity_ge{cutoff_label.replace('.', 'p')}")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix, table, include_index in [
        ("matrix", matrix, True),
        ("pair_counts", counts, True),
        ("genome_pairs", pairs, False),
    ]:
        path = Path(f"{prefix}_{suffix}.tsv")
        table.to_csv(path, sep="\t", index=include_index, na_rep="NA")
        outputs.append(path)
    fig = make_matrix_figure(
        matrix, summary,
        title=f"Mean RBM fraction at >= {cutoff_label}% sequence identity",
        colorbar_label=f"Mean RBM fraction (identity >= {cutoff_label}%)",
        vmin=0, vmax=1,
        footnotes=(
            "Cells: mean per-genome-pair ratio; outlined diagonal: within-genomovar comparisons.",
            f"Ratio = RBM rows with column 3 >= {cutoff_label}% / all RBM rows; "
            "missing comparisons excluded.",
        ),
    )
    for extension in ["png", "pdf", "svg"]:
        path = Path(f"{prefix}.{extension}")
        fig.savefig(path, dpi=450, facecolor="white")
        outputs.append(path)
    plt.close(fig)
    available = pairs.loc[pairs["status"] == "available"]
    report_path = Path(f"{prefix}_report.txt")
    outputs.append(report_path)
    report_lines = [
        "RBM sequence-identity fraction matrix",
        f"RBM directory: {rbm_dir.resolve()}",
        f"Inclusive identity cutoff: {cutoff_label}%",
        f"Seed used by the shared ANI direction selection: {seed}",
        f"RBM files found: {inventory['files_found']}",
        f"Genome self-comparison files excluded: {inventory['self_files_excluded']}",
        f"Out-of-scope files excluded: {inventory['out_of_scope_files_excluded']}",
        f"Eligible reciprocal file pairs reduced to one: {inventory['reciprocal_pairs']}",
        f"Genomovars: {len(summary)}",
        f"Genomes: {int(summary['n_genomes'].sum())}",
        f"Expected unique genome pairs: {len(pairs)}",
        f"Available unique genome pairs: {len(available)}",
        f"Missing unique genome pairs: {len(pairs) - len(available)}",
        f"Total RBM rows across selected files: {int(available['total_rows'].sum())}",
        f"Qualifying rows across selected files: {int(available['qualifying_rows'].sum())}",
        "Axis order: " + ", ".join(summary["genomovar_complete"]),
        "Phylogroup order: " + ", ".join(summary["phylogroup"].drop_duplicates()),
        "",
        "Per-file ratio: count(column 3 >= cutoff) / count(all nonblank RBM rows).",
        "Each reported row is counted once; no length weighting or gene-row deduplication.",
        "Denominator is the number of reported RBMs, NOT the total genes in either genome.",
        "Identity is interpreted at the precision reported in column 3; reported 100% "
        "does not establish identity across full-length genes.",
        "Genomovar cells: arithmetic mean of per-genome-pair ratios, with equal "
        "weight per genome pair, NOT pooled qualifying counts / pooled total counts.",
        "Diagonal: mean ratios for distinct genomes in the same genomovar.",
        "Reciprocal files: use the retained ANI direction if both exist; use the "
        "available direction if only one exists. Unused files are recorded in the pair table.",
        "Missing pairs are reported and omitted from means, never assigned ratio zero.",
        "An empty matrix cell is NA; pair-count matrices count available genome pairs.",
        "Axes, phylogroup annotations, and membership match the ANI heatmap.",
        "Heatmap scale is 0 to 1; annotations use 2 decimals; TSV ratios are unrounded.",
        f"Input paths, ANI exclusions, and software versions: {output_prefix}_report.txt",
        "",
        *[f"WARNING: {warning}" for warning in warnings],
        "Outputs:",
        *[str(path.resolve()) for path in outputs],
    ]
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("\n" + "\n".join(report_lines))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ani", type=Path,
                        default=SCRIPT_DIR.parent / "inputs" / "EcoZUR.ani",
                        help="Five-column fastANI input.")
    parser.add_argument("--metadata", type=Path,
                        default=SCRIPT_DIR.parent / "inputs" / "EcoZUR_genomovar_metadata.tsv",
                        help="TSV with genome, genomovar_complete, and phylogroup.")
    parser.add_argument("--output-prefix", type=Path,
                        default=SCRIPT_DIR.parent / "results" / "within_genomovar_ANI",
                        help="Output path prefix (no extension; existing outputs overwritten).")
    parser.add_argument("--seed", type=nonnegative_int, default=13,
                        help="Random seed for tied fragment counts and point jitter (default: 13).")
    parser.add_argument("--rbm", action="store_true",
                        help="Also generate the RBM identity-ratio heatmap (default cutoff: 100%%).")
    parser.add_argument("--rbm-dir", type=Path,
                        default=SCRIPT_DIR.parent / "inputs" / "RBM",
                        help="Directory of genomeA___genomeB.rbm files; used only when RBM analysis is enabled.")
    parser.add_argument("--rbm-min-identity", type=identity_percentage,
                        help="Enable RBM analysis at this inclusive percent-identity cutoff, e.g. 99.8.")
    args = parser.parse_args(argv)
    ani = read_ani(args.ani)
    metadata = read_metadata(args.metadata)
    input_genomes = set(ani["query"]) | set(ani["reference"])
    selected, tied_pairs = select_pairs(ani, args.seed)
    within, summary, warnings = summarize_groups(selected, metadata, input_genomes)
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    fig = make_figure(within, summary, args.seed)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix, table in [
        ("selected_pairs.tsv", selected),
        ("within_pairs.tsv", within),
        ("summary.tsv", summary),
    ]:
        output = Path(f"{args.output_prefix}_{suffix}")
        table.to_csv(output, sep="\t", index=False)
        outputs.append(output)
    for extension in ["png", "pdf", "svg"]:
        output = Path(f"{args.output_prefix}.{extension}")
        fig.savefig(output, dpi=450, facecolor="white")
        outputs.append(output)
    plt.close(fig)
    report_lines = [
        f"ANI input: {args.ani.resolve()}",
        f"Metadata: {args.metadata.resolve()}",
        f"Seed: {args.seed}",
        f"Input ANI rows: {len(ani)}",
        f"Unique input genomes: {len(input_genomes)}",
        f"Self comparisons removed: {int(ani['query'].eq(ani['reference']).sum())}",
        f"Selected unordered pairs (before metadata filtering): {len(selected)}",
        f"Equal-maximum fragment-count pairs randomized: {tied_pairs}",
        f"Metadata genomes: {len(metadata)}",
        f"Genomovars plotted (at least {MIN_MEMBERS} input genomes): {len(summary)}",
        f"Genomes in plotted genomovars: {int(summary['n_genomes'].sum())}",
        f"Within-genomovar pairs plotted: {len(within)}",
        "Genomovar order: " + ", ".join(summary["genomovar_complete"]),
        "Phylogroup order: " + ", ".join(summary["phylogroup"].drop_duplicates()),
        "",
        "Method: retain the entire row with maximum column-5 query fragments per "
        "unordered pair; select uniformly at random among tied maximum rows.",
        "ANI summaries: unweighted descriptive statistics of retained within-group "
        "pairwise ANI percentages; sample standard deviation (ddof=1).",
        "Membership: distinct annotated genomes present in the ANI input.",
        "No reciprocal ANI averaging; no ANI values imputed or rounded before analysis.",
        "Violin: Gaussian KDE, Scott bandwidth, cut=0, equal maximum widths; "
        "median bar; seeded uniform horizontal point jitter.",
        "Pairwise ANI observations share genomes and are not independent replicates. "
        "Violin density estimates are especially uncertain for small genomovars.",
        "",
        *[f"WARNING: {warning}" for warning in warnings],
        "",
        f"Versions: numpy={np.__version__}, pandas={pd.__version__}, "
        f"matplotlib={matplotlib.__version__}, seaborn={sns.__version__}",
    ]
    report = Path(f"{args.output_prefix}_report.txt")
    outputs.append(report)
    report_lines.extend(["", "Outputs:", *[str(path.resolve()) for path in outputs]])
    report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("\n".join(report_lines))
    write_genomovar_matrix(selected, metadata, summary, args.output_prefix, args.seed)
    if args.rbm or args.rbm_min_identity is not None:
        threshold = 100.0 if args.rbm_min_identity is None else args.rbm_min_identity
        write_rbm_matrix(selected, metadata, summary, input_genomes,
                         args.rbm_dir, args.output_prefix, threshold, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
