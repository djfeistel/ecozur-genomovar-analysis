#!/usr/bin/env python3
"""Compare reciprocal fastANI values (A->B vs B->A) for every genome pair."""

import argparse
import math
import random
import statistics
import sys

DEFAULT_INPUT = "../input/EcoZUR.ani"

EPILOG = """
what it does:
  fastANI is not symmetric: the ANI of query A vs reference B can differ
  slightly from query B vs reference A. For every unordered pair {A, B}
  found in BOTH directions, this script computes
      |ANI(A -> B) - ANI(B -> A)|
  and prints the mean and standard deviation of those differences
  (plus n, median, min and max). Self comparisons (A vs A) are ignored.

  Example:
      B001-5  E21-3   98.7222  1435  1665
      E21-3   B001-5  98.7676  1438  1613
      -> |98.7222 - 98.7676| = 0.0454

input format (fastANI output, tab-separated, no header):
  query  reference  ANI  matched_fragments  total_fragments

plot (-p):
  (A) violin + boxplot of all |differences| (red diamond = mean), showing
      how large the asymmetry typically is and how heavy the tail is.
  (B) violins of |difference| split by the pair's mean ANI
      ((ANI(A->B) + ANI(B->A)) / 2), binned at --ani-bins, showing whether
      the asymmetry depends on how similar the genomes are.
  --points overlays every pair as a grey dot (alpha 0.1), jittered
  horizontally in proportion to the violin width, so the raw spread is visible.

used-vs-unused plot (-u):
  Downstream methods keep, for each pair, the row with the LARGER fifth
  column (total query fragments, i.e. the larger genome as query); ties are
  broken uniformly at random with --seed. This scatter plots the discarded
  direction (smaller genome as query; x) against the retained one (y), with
  the y = x diagonal and an ordinary least-squares trend line (slope,
  intercept and R^2 in the legend). Requires the fifth column.

  Output format follows the file extension (.png, .pdf, .svg). Plots require
  matplotlib.

examples:
  python ani_reciprocal_difference.py
  python ani_reciprocal_difference.py -i other_run.ani
  python ani_reciprocal_difference.py -o reciprocal_differences.tsv
  python ani_reciprocal_difference.py -p ../results/reciprocal_differences.png
  python ani_reciprocal_difference.py -p diffs.pdf --ani-bins 97,98,99,99.5
  python ani_reciprocal_difference.py -p diffs.png --points
  python ani_reciprocal_difference.py -u ../results/used_vs_unused_ani.png
"""


def parse_args():
    p = argparse.ArgumentParser(
        prog="ani_reciprocal_difference.py",
        description=("For each genome pair, compute the absolute difference between "
                     "the fastANI value in one direction (A vs B) and the reverse "
                     "(B vs A), then print the mean and SD over all pairs."),
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", default=DEFAULT_INPUT,
                   help="fastANI output file (default: %(default)s)")
    p.add_argument("-o", "--output", metavar="TSV",
                   help="Optional: also write every pair and its difference to this TSV "
                        "(columns: genome_A, genome_B, ANI_A_vs_B, ANI_B_vs_A, abs_diff)")
    p.add_argument("-p", "--plot", metavar="FILE",
                   help="Optional: save a violin-plot figure of the differences to FILE "
                        "(.png, .pdf, .svg, ...)")
    p.add_argument("--ani-bins", default="97,98,99,99.5", metavar="LIST",
                   help="Comma-separated mean-ANI breakpoints for plot panel B "
                        "(default: %(default)s)")
    p.add_argument("--points", action="store_true",
                   help="Overlay every pair as a grey dot (alpha 0.1) on the -p violins")
    p.add_argument("-u", "--used-plot", metavar="FILE",
                   help="Optional: save a scatter of unused (x) vs used (y) ANI per pair, "
                        "with y = x diagonal and linear regression trend line, to FILE")
    p.add_argument("--seed", type=int, default=13,
                   help="Random seed for fragment-count ties and point jitter "
                        "(default: %(default)s)")
    args = p.parse_args()
    try:
        args.ani_bins = sorted(float(x) for x in args.ani_bins.split(",") if x.strip())
    except ValueError:
        p.error(f"--ani-bins must be comma-separated numbers, got {args.ani_bins!r}")
    if args.points and not args.plot:
        p.error("--points only applies to the violin plot; also give -p FILE")
    return args


def read_ani(path):
    """Return {(query, ref): ANI}, {(query, ref): total query fragments or None}, n duplicates."""
    ani, frags = {}, {}
    duplicates = 0
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                sys.exit(f"Error: line {n} has fewer than 3 tab-separated columns: {line!r}")
            query, ref = fields[0], fields[1]
            try:
                value = float(fields[2])
            except ValueError:
                sys.exit(f"Error: line {n}: ANI value is not a number: {fields[2]!r}")
            total = None
            if len(fields) >= 5:
                try:
                    total = int(fields[4])
                except ValueError:
                    sys.exit(f"Error: line {n}: total_fragments is not an integer: {fields[4]!r}")
            if (query, ref) in ani:
                duplicates += 1
            ani[(query, ref)] = value
            frags[(query, ref)] = total
    return ani, frags, duplicates


def split_used_unused(pairs, frags, seed):
    """Per pair, the retained (larger query fragments; seeded random tie) and discarded ANI."""
    rng = random.Random(seed)
    used, unused, ties = [], [], 0
    for a, b, ab, ba, _ in pairs:
        fa, fb = frags[(a, b)], frags[(b, a)]
        if fa is None or fb is None:
            sys.exit("Error: -u needs the fifth column (total query fragments) for every row")
        if fa == fb:
            ties += 1
            keep_ab = rng.random() < 0.5
        else:
            keep_ab = fa > fb
        used.append(ab if keep_ab else ba)
        unused.append(ba if keep_ab else ab)
    return used, unused, ties


def load_pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Error: plotting requires matplotlib (pip/conda install matplotlib)")
    return plt


def bin_label(lo, hi):
    if lo is None:
        return f"<{hi:g}"
    if hi is None:
        return f"\u2265{lo:g}"
    return f"{lo:g}\u2013{hi:g}"


def _kde(values, n_grid=200):
    """Gaussian KDE (Scott's bandwidth, as matplotlib's violinplot) on a grid
    spanning the data; vectorized in chunks because matplotlib's own KDE is
    very slow for tens of thousands of points."""
    import numpy as np
    v = np.asarray(values, dtype=float)
    grid = np.linspace(v.min(), v.max(), n_grid)
    bw = v.std(ddof=1) * v.size ** (-1 / 5)
    dens = np.zeros(n_grid)
    for i in range(0, v.size, 5000):
        z = (grid[:, None] - v[None, i:i + 5000]) / bw
        dens += np.exp(-0.5 * z * z).sum(axis=1)
    return grid, dens / (v.size * bw * math.sqrt(2 * math.pi))


def _violin_box(ax, groups, positions, rng=None, width=0.8):
    """Violin (where a KDE is possible) + narrow boxplot + mean marker for each group.

    If rng is given, every value is also drawn as a grey dot, jittered
    horizontally in proportion to the violin's width at that value.
    """
    import numpy as np
    stats, stat_pos = [], []
    for g, x in zip(groups, positions):
        kde_ok = len(g) > 1 and min(g) != max(g)
        if kde_ok:
            grid, dens = _kde(g)
            stats.append({"coords": grid, "vals": dens, "mean": statistics.fmean(g),
                          "median": statistics.median(g), "min": min(g), "max": max(g)})
            stat_pos.append(x)
        if rng is not None:
            half = width / 2
            scale = np.interp(g, grid, dens / dens.max()) if kde_ok else 0.3
            xs = x + rng.uniform(-1, 1, len(g)) * half * scale
            ax.scatter(xs, g, s=3, color="grey", alpha=0.1, linewidths=0,
                       zorder=2, rasterized=True)
    if stats:
        parts = ax.violin(stats, positions=stat_pos, widths=width, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor("#7aa6c2")
            body.set_edgecolor("#2f5d7c")
            body.set_alpha(0.7)
            body.set_zorder(1)
    ax.boxplot(groups, positions=positions, widths=0.12, showfliers=False,
               patch_artist=True, medianprops={"color": "black"},
               boxprops={"facecolor": "white", "edgecolor": "black"}, zorder=3)
    ax.scatter(positions, [statistics.fmean(g) for g in groups], marker="D",
               color="#c0392b", s=18, zorder=4, label="mean")
    # Mean/SD label just above y = 0, right-aligned to end just left of the box
    # (x in data units, y in axes units so it stays at the bottom for any ylim).
    for g, x in zip(groups, positions):
        sd = statistics.stdev(g) if len(g) > 1 else math.nan
        ax.text(x - 0.09, 0.015, f"mean {statistics.fmean(g):.4f}\nSD {sd:.4f}",
                transform=ax.get_xaxis_transform(), ha="right", va="bottom",
                fontsize=7, linespacing=1.15, zorder=5,
                bbox={"boxstyle": "round,pad=0.2", "facecolor": "white",
                      "edgecolor": "none", "alpha": 0.75})


def make_plot(pairs, path, breakpoints, points=False, seed=13):
    plt = load_pyplot()
    rng = None
    if points:
        import numpy as np
        rng = np.random.default_rng(seed)

    diffs = [p[4] for p in pairs]
    edges = [None] + breakpoints + [None]
    labels, groups = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        g = [p[4] for p in pairs
             if (lo is None or (p[2] + p[3]) / 2 >= lo) and (hi is None or (p[2] + p[3]) / 2 < hi)]
        if g:
            labels.append(bin_label(lo, hi))
            groups.append(g)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True,
                                   gridspec_kw={"width_ratios": [1, max(len(groups), 2)]})

    _violin_box(ax1, [diffs], [1], rng)
    ax1.set_xticks([1], [f"all pairs\n(n={len(diffs):,})"])
    ax1.set_ylabel("|ANI(A vs B) \u2212 ANI(B vs A)|  (ANI %)")
    ax1.set_title("A. All pairs", loc="left", fontsize=11)
    ax1.legend(loc="upper right", frameon=False, fontsize=9)

    pos = list(range(1, len(groups) + 1))
    _violin_box(ax2, groups, pos, rng)
    ax2.set_xticks(pos, [f"{lab}\n(n={len(g):,})" for lab, g in zip(labels, groups)])
    ax2.set_xlabel("Mean reciprocal ANI of the pair (%)")
    ax2.set_title("B. By ANI level", loc="left", fontsize=11)

    for ax in (ax1, ax2):
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0.42)
        ax.grid(axis="y", alpha=0.3)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def make_used_plot(used, unused, fit, path):
    plt = load_pyplot()
    slope, intercept, r2 = fit
    sign = "+" if intercept >= 0 else "\u2212"

    lo = min(min(used), min(unused))
    hi = max(max(used), max(unused))
    pad = (hi - lo) * 0.02 or 0.1
    lim = (lo - pad, hi + pad)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(unused, used, s=4, color="grey", alpha=0.2, linewidths=0,
               rasterized=True, zorder=1, label=f"genome pairs (n={len(used):,})")
    ax.plot(lim, lim, ls="--", color="black", lw=1, zorder=2, label="y = x")
    ax.plot(lim, [slope * x + intercept for x in lim], color="#c0392b", lw=1.5, zorder=3,
            label=(f"linear fit: y = {slope:.4f}x {sign} "
                   f"{abs(intercept):.4f}\nR\u00b2 = {r2:.4f}"))
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Unused ANI (%)  [smaller genome as query]")
    ax.set_ylabel("Used ANI (%)  [larger genome as query]")
    ax.set_title("Retained vs discarded fastANI direction", fontsize=11)
    ax.grid(alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", frameon=False, fontsize=9, markerscale=3)

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main():
    args = parse_args()
    try:
        ani, frags, duplicates = read_ani(args.input)
    except FileNotFoundError:
        sys.exit(f"Error: input file not found: {args.input}")

    pairs = []
    one_way_only = 0
    for (a, b), ab in ani.items():
        if a == b:
            continue
        ba = ani.get((b, a))
        if ba is None:
            one_way_only += 1
        elif a < b:
            pairs.append((a, b, ab, ba, abs(ab - ba)))

    if not pairs:
        sys.exit("Error: no reciprocal (A vs B and B vs A) comparisons found")

    diffs = [p[4] for p in pairs]
    mean = statistics.fmean(diffs)
    sd = statistics.stdev(diffs) if len(diffs) > 1 else math.nan
    genomes = {g for key in ani for g in key}
    worst = max(pairs, key=lambda p: p[4])

    print(f"Input file:                     {args.input}")
    print(f"Genomes:                        {len(genomes)}")
    print(f"Comparisons read (rows):        {len(ani)}")
    print(f"Reciprocal pairs compared:      {len(pairs)}")
    if one_way_only:
        print(f"Comparisons missing a reverse:  {one_way_only} (skipped)")
    if duplicates:
        print(f"Duplicate rows (last kept):     {duplicates}")
    print()
    print("|ANI(A vs B) - ANI(B vs A)|")
    print(f"  Mean:    {mean:.4f}")
    print(f"  SD:      {sd:.4f}")
    print(f"  Median:  {statistics.median(diffs):.4f}")
    print(f"  Min:     {min(diffs):.4f}")
    print(f"  Max:     {worst[4]:.4f}  ({worst[0]} vs {worst[1]}: {worst[2]} / {worst[3]})")

    if args.output:
        with open(args.output, "w") as fh:
            fh.write("genome_A\tgenome_B\tANI_A_vs_B\tANI_B_vs_A\tabs_diff\n")
            for a, b, ab, ba, d in sorted(pairs, key=lambda p: -p[4]):
                fh.write(f"{a}\t{b}\t{ab}\t{ba}\t{d:.4f}\n")
        print(f"\nPer-pair differences written to: {args.output}")

    if args.plot:
        make_plot(pairs, args.plot, args.ani_bins, args.points, args.seed)
        print(f"\nPlot written to: {args.plot}")

    if args.used_plot:
        used, unused, ties = split_used_unused(pairs, frags, args.seed)
        if len(used) < 2 or len(set(unused)) < 2:
            sys.exit("Error: -u needs at least two pairs with differing unused ANI for a regression")
        slope, intercept = statistics.linear_regression(unused, used)
        r2 = statistics.correlation(unused, used) ** 2
        signed = [u - x for u, x in zip(used, unused)]
        print()
        print("Used (larger genome as query) vs unused (smaller genome as query) ANI")
        print(f"  Fragment-count ties:     {ties} (broken at random, seed {args.seed})")
        print(f"  Mean used - unused:      {statistics.fmean(signed):+.4f}")
        print(f"  Used > unused:           {sum(s > 0 for s in signed)} of {len(signed)}")
        print(f"  Linear fit (y=used):     y = {slope:.4f}x {intercept:+.4f},  R^2 = {r2:.4f}")
        make_used_plot(used, unused, (slope, intercept, r2), args.used_plot)
        print(f"\nUsed-vs-unused plot written to: {args.used_plot}")


if __name__ == "__main__":
    main()
