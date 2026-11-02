#!/usr/bin/env python3
"""Scatter-plot any two statistics from an official CheckM bin_stats file
(bin_stats_ext.tsv, bin_stats.analyze.tsv or bin_stats.tree.tsv)."""

import argparse
import ast
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path
from matplotlib.ticker import AutoMinorLocator, Formatter, NullFormatter

DEFAULT_INPUT = "bin_stats_ext.tsv"
NAME_COL = "genome"
# Okabe-Ito colour-blind-safe palette for points; red tint for failing regions
POINT_COLOR = "#0072B2"
HIGHLIGHT_COLOR = "#E69F00"
SHADE_COLOR = "#E31A1C"
# ggplot2 theme_grey colours
GREY20, GREY30, GREY92 = "#333333", "#4D4D4D", "#EBEBEB"
# Liberation Sans is metrically identical to Arial (standard for journals)
FONTS = ["Arial", "Liberation Sans", "Helvetica", "Nimbus Sans", "DejaVu Sans"]
# Publication-style axis titles: column -> (label, unit). Columns in BP_COLUMNS
# get bp/kbp/Mbp chosen automatically from the data range.
PRETTY = {
    "# genomes": ("No. of reference genomes", None),
    "# markers": ("No. of marker genes", None),
    "# marker sets": ("No. of marker sets", None),
    "0 copies": ("Markers with 0 copies", None),
    "1 copy": ("Markers with 1 copy", None),
    "2 copies": ("Markers with 2 copies", None),
    "3 copies": ("Markers with 3 copies", None),
    "4 copies": ("Markers with 4 copies", None),
    "5+ copies": ("Markers with 5+ copies", None),
    "Completeness": ("Completeness", "%"),
    "Contamination": ("Contamination", "%"),
    "GC": ("GC content", None),
    "GC std": ("GC content SD", None),
    "Genome size": ("Genome size", "bp"),
    "# ambiguous bases": ("No. of ambiguous bases", None),
    "# scaffolds": ("No. of scaffolds", None),
    "# contigs": ("No. of contigs", None),
    "Longest scaffold": ("Longest scaffold", "bp"),
    "Longest contig": ("Longest contig", "bp"),
    "N50 (scaffolds)": ("Scaffold N50", "bp"),
    "N50 (contigs)": ("Contig N50", "bp"),
    "Mean scaffold length": ("Mean scaffold length", "bp"),
    "Mean contig length": ("Mean contig length", "bp"),
    "Coding density": ("Coding density", None),
    "Translation table": ("Translation table", None),
    "# predicted genes": ("No. of predicted genes", None),
}
# CheckM names the marker copy-number counts "0".."5+"; bare digits would clash
# with selecting columns by number, so give them descriptive names.
RENAME = {"0": "0 copies", "1": "1 copy", "2": "2 copies",
          "3": "3 copies", "4": "4 copies", "5+": "5+ copies"}

EPILOG = """
input file format:
  Official CheckM (v1) per-bin statistics, one bin per line:
      <bin name><TAB>{'Completeness': 100.0, 'Contamination': 0.48, ...}
  Accepted files (all found in <checkm_output>/storage/):
      bin_stats_ext.tsv      (DEFAULT) everything: quality + assembly stats
      bin_stats.analyze.tsv  assembly stats only (NO Completeness/Contamination)
      bin_stats.tree.tsv     assembly stats only (NO Completeness/Contamination)
  The per-marker ID lists (GCN0 ... GCN5+) are not plottable and are skipped;
  their counts are available as the "0 copies" ... "5+ copies" columns.

available columns in bin_stats_ext.tsv (name or number may be used):
   1  genome               (bin name - used for --label / --exclude)
   2  marker lineage       (text, e.g. f__Enterobacteriaceae - not plottable)
   3  # genomes            (reference genomes behind the marker set)
   4  # markers            (marker genes in the lineage-specific set)
   5  # marker sets        (collocated marker sets)
   6  0 copies             (markers found 0 times = missing)
   7  1 copy               (markers found once)
   8  2 copies             (markers found twice)
   9  3 copies
  10  4 copies
  11  5+ copies
  12  Completeness         (%)
  13  Contamination        (%)
  14  GC                   (fraction 0-1)
  15  GC std
  16  Genome size          (bp)
  17  # ambiguous bases
  18  # scaffolds
  19  # contigs
  20  Longest scaffold     (bp)
  21  Longest contig       (bp)
  22  N50 (scaffolds)      (bp)
  23  N50 (contigs)        (bp)
  24  Mean scaffold length (bp)
  25  Mean contig length   (bp)
  26  Coding density       (fraction 0-1)
  27  Translation table
  28  # predicted genes

  bin_stats.analyze.tsv / bin_stats.tree.tsv contain only columns 14-28 (in
  that order), so their column NUMBERS differ: ALWAYS run --list first when
  using a different -i file, or select columns by name.

  Column names are matched case-insensitively. Names containing spaces or '#'
  MUST be quoted, e.g. -y "Genome size"  or  -x "# scaffolds".
  A bare number always means a column NUMBER (so -x 7 = "1 copy").

output:
  Always a PNG at 300 dpi. If -o is not given the file is named
      CheckM_<Y>_vs_<X>.png
  (spaces/special characters replaced by '_'), e.g.
      CheckM_Contamination_vs_Completeness.png
  If -o is given without a .png extension, .png is appended.

examples:
  # Classic CheckM quality plot
  python plot_checkm_stats.py -x Completeness -y Contamination

  # Use column numbers instead of names (16 = Genome size, 14 = GC)
  python plot_checkm_stats.py -x 16 -y 14

  # Genome size vs N50, log-scaled y, with genome names on the points
  python plot_checkm_stats.py -x "Genome size" -y "N50 (scaffolds)" --logy --label

  # Use the assembly-only analyze file instead
  python plot_checkm_stats.py -i bin_stats.analyze.tsv -x "Genome size" \\
        -y "Coding density" --label

  # Colour points by contamination and save under a custom name
  python plot_checkm_stats.py -x Completeness -y "# scaffolds" \\
        --color-by Contamination -o my_plot.png

  # Zoom in: X axis from 84 to 100, Y axis from 0 to 10
  python plot_checkm_stats.py -x Completeness -y Contamination --xlim 84 100 --ylim 0 10

  # Shade the failing region red: Completeness < 90 OR Contamination > 5
  # (MIMAG high-quality thresholds). Genomes inside the red tint fail.
  python plot_checkm_stats.py -x Completeness -y Contamination --xline 90 --yline 5

  # Two quality tiers: the region failing even 50/10 gets a darker red
  python plot_checkm_stats.py -x Completeness -y Contamination --xline 50 90 --yline 5 10

  # Publication figure: no title, custom axis labels, larger fonts
  python plot_checkm_stats.py -x Completeness -y Contamination --xline 90 --yline 5 \\
        --no-title --base-size 13 --xlabel "CheckM completeness (%)"

plot style:
  The figure uses the default ggplot2 look (theme_grey): light grey panel with
  no border, solid white major grid lines, dashed white minor grid lines,
  outward ticks, Arial-compatible font. Axis labels get units
  automatically, e.g. "Completeness (%)", "Genome size (Mbp)", "Scaffold N50
  (kbp)". Base-pair columns are DISPLAYED in kbp/Mbp, but --xlim, --ylim,
  --xline and --yline always use the raw values in the file (bp, %, 0-1 GC).
  --base-size PT   base font size (default 11; axis ticks are 80% of it)
  --xlabel/--ylabel TEXT   override the automatic axis titles
  --no-title       omit the plot title (common for journal figures)
  --no-legend      omit the legend (normally centred under the plot)

axis limits:
  --xlim MIN MAX   X axis range, e.g. --xlim 84 100
  --ylim MIN MAX   Y axis range, e.g. --ylim 0 10
  --pad FRACTION   margin added to both ends of --xlim/--ylim so dots on the
                   limits are fully visible (default 0.02 = 2% of the range;
                   --pad 0 gives exact limits, --pad 0.05 gives more room)

threshold shading (--xline / --yline):
  --xline VALUE    X threshold, e.g. --xline 90 for Completeness
  --yline VALUE    Y threshold, e.g. --yline 5 for Contamination
  The region where a genome FAILS the thresholds is tinted red (alpha 0.2);
  points inside the tint are outside the accepted range. Which side fails is
  set by --x-pass / --y-pass (see below). With the defaults this is
      X < xline  (e.g. Completeness < 90)   OR   Y > yline  (e.g. Contamination > 5)
  i.e. an L-shaped region left of the X threshold and above the Y threshold;
  the white rectangle that remains is the passing area. Giving only --xline
  or only --yline shades just that side. If several values are given, each
  extra (more lenient) tier is drawn on top, so worse genomes sit in darker red.
  The axes are automatically widened to include the thresholds.
  --shade-color COLOR   tint colour (default #E31A1C, red)
  --shade-alpha A       tint transparency 0-1 (default 0.2)
  --shade-label TEXT    legend text for the tint (default: the fail rule,
                        e.g. "Completeness < 90 or Contamination > 5")
  --label-failed        write the genome NAME next to every point that falls
                        in the red (failing) region. Labels are placed like
                        ggrepel: each one is moved to the nearest free spot
                        so it overlaps no data point and no other label, stays
                        inside the plot, and gets a thin leader line back to
                        its point when moved away. Genomes outside a zoomed
                        --xlim/--ylim range are not labelled. Combined with
                        --exclude FILE --highlight-excluded, the highlighted
                        (orange diamond) genomes are labelled too, so it also
                        works without --xline/--yline.
      python plot_checkm_stats.py -x Completeness -y Contamination \\
            --xline 90 --yline 5 --label-failed
      python plot_checkm_stats.py -x "Genome size" -y "# scaffolds" \\
            --exclude CheckM_failed_Completeness_ge90_Contamination_le5.txt \\
            --highlight-excluded --label-failed

filtering genomes with the --xline / --yline thresholds:
  --write-failed [FILE]
      Treat --xline/--yline as MINIMUM QUALITY CRITERIA and write the names of
      every genome that FAILS them to a text file (one name per line). The
      plot is still made as usual. If FILE is omitted the file is named
          CheckM_failed_<X>_<rule><value>_<Y>_<rule><value>.txt
      e.g. CheckM_failed_Completeness_ge90_Contamination_le5.txt
  --x-pass ge|le   a genome passes the X criterion if X >= (ge) or <= (le)
                   the --xline value (default: ge, e.g. Completeness >= 90)
  --y-pass ge|le   same for Y (default: le, e.g. Contamination <= 5)
  These rules also decide which side of each threshold is shaded red.
  A genome FAILS if it misses ANY criterion, or if its value is missing.
  Only the axes that have a line are used (e.g. only --xline -> X only).
  If several values are given to --xline/--yline, the strictest one is used.

  --exclude FILE   Remove the genomes listed in FILE from the data BEFORE
                   plotting (or before --write-failed). FILE = one genome name
                   per line (first word on the line is used; blank lines and
                   lines starting with '#' are ignored). Can be repeated to
                   combine several lists.
  --highlight-excluded
                   With --exclude, keep listed genomes and draw them as orange
                   diamonds instead of removing them. Other genomes remain
                   circles; a legend identifies both groups. Requires --exclude.
                   With --color-by, only the circles use the colour scale.
                   Missing plotting values are still skipped as usual.
                   With --write-failed, retained highlighted genomes are also
                   evaluated against the thresholds.

  typical workflow:
    # 1) QC plot: write genomes with Completeness < 90 or Contamination > 5
    python plot_checkm_stats.py -x Completeness -y Contamination \\
          --xline 90 --yline 5 --write-failed
    # 2) replot any other columns without those genomes
    python plot_checkm_stats.py -x "Genome size" -y "# scaffolds" \\
          --exclude CheckM_failed_Completeness_ge90_Contamination_le5.txt
    # 3) show the listed genomes instead, using orange diamonds
    python plot_checkm_stats.py -x "Genome size" -y "# scaffolds" \\
          --exclude CheckM_failed_Completeness_ge90_Contamination_le5.txt \\
          --highlight-excluded -o highlighted.png
"""

AUTO = "__auto__"


def parse_args():
    p = argparse.ArgumentParser(
        prog="plot_checkm_stats.py",
        description=(
            "Make a scatter plot of any two statistics from an official CheckM "
            "bin_stats file (default: bin_stats_ext.tsv). Each point is one genome. "
            "Choose the X and Y columns by name or by 1-based column number "
            "(see --list)."
        ),
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", default=DEFAULT_INPUT,
                   help="CheckM bin_stats file: bin_stats_ext.tsv, bin_stats.analyze.tsv "
                        "or bin_stats.tree.tsv (default: %(default)s)")
    p.add_argument("-x", "--x", dest="xcol",
                   help="Column for the X axis (header name or column number)")
    p.add_argument("-y", "--y", dest="ycol",
                   help="Column for the Y axis (header name or column number)")
    p.add_argument("-o", "--output",
                   help="Output PNG file (default: CheckM_<Y>_vs_<X>.png)")
    p.add_argument("--list", action="store_true",
                   help="Print the available columns of the input file and exit")
    p.add_argument("--color-by", metavar="COLUMN",
                   help="Optional numeric column used to colour the points (adds a colour bar)")
    p.add_argument("--label", action="store_true",
                   help="Write the genome name next to every point")
    p.add_argument("--label-failed", action="store_true",
                   help="Label the genomes that fail --xline/--yline (the points in the "
                        "red region) and, with --highlight-excluded, the highlighted "
                        "genomes, using non-overlapping labels and leader lines")
    p.add_argument("--logx", action="store_true", help="Log-scale the X axis")
    p.add_argument("--logy", action="store_true", help="Log-scale the Y axis")
    p.add_argument("--xlim", type=float, nargs=2, metavar=("MIN", "MAX"),
                   help="Set the X axis range, e.g. --xlim 84 100 (default: automatic)")
    p.add_argument("--ylim", type=float, nargs=2, metavar=("MIN", "MAX"),
                   help="Set the Y axis range, e.g. --ylim 0 10 (default: automatic)")
    p.add_argument("--pad", type=float, default=0.02, metavar="FRACTION",
                   help="Extra margin added to each end of --xlim/--ylim as a fraction "
                        "of the range, so points sitting exactly on MIN or MAX are not "
                        "cut off (default: %(default)s = 2%%; use 0 for exact limits)")
    p.add_argument("--xline", "--vline", dest="xline", type=float, nargs="+",
                   action="extend", metavar="X",
                   help="X threshold: tint the failing side (default: X < value) red "
                        "(one or more values, e.g. --xline 90 or --xline 50 90)")
    p.add_argument("--yline", "--hline", dest="yline", type=float, nargs="+",
                   action="extend", metavar="Y",
                   help="Y threshold: tint the failing side (default: Y > value) red "
                        "(one or more values, e.g. --yline 5 or --yline 5 10)")
    p.add_argument("--shade-color", default=SHADE_COLOR, metavar="COLOR",
                   help="Colour of the failing-region tint (default: %(default)s)")
    p.add_argument("--shade-alpha", type=float, default=0.2, metavar="A",
                   help="Transparency of the failing-region tint, 0-1 (default: %(default)s)")
    p.add_argument("--shade-label", metavar="TEXT",
                   help="Legend text for the tint (default: the fail rule, e.g. "
                        "'Completeness < 90 or Contamination > 5')")
    p.add_argument("--write-failed", nargs="?", const=AUTO, metavar="FILE",
                   help="Use --xline/--yline as minimum criteria and write the names of "
                        "genomes that FAIL them to FILE (default name: "
                        "CheckM_failed_<X>_<rule><value>_<Y>_<rule><value>.txt)")
    p.add_argument("--x-pass", choices=["ge", "le"], default="ge",
                   help="Pass rule for X (shading and --write-failed): ge = X >= --xline, "
                        "le = X <= --xline (default: %(default)s)")
    p.add_argument("--y-pass", choices=["ge", "le"], default="le",
                   help="Pass rule for Y (shading and --write-failed): ge = Y >= --yline, "
                        "le = Y <= --yline (default: %(default)s)")
    p.add_argument("--exclude", action="append", metavar="FILE",
                   help="Remove the genomes listed in FILE (one name per line, e.g. a "
                        "--write-failed file) before plotting; can be repeated. "
                        "Use --highlight-excluded to highlight them instead")
    p.add_argument("--highlight-excluded", action="store_true",
                   help="Keep genomes listed by --exclude and plot them as orange "
                        "diamonds with a legend instead of removing them "
                        "(requires --exclude)")
    p.add_argument("--title", help="Plot title (default: '<Y> vs <X> (n=<genomes>)')")
    p.add_argument("--no-title", action="store_true", help="Do not draw a plot title")
    p.add_argument("--xlabel", metavar="TEXT",
                   help="X axis title (default: column name with units)")
    p.add_argument("--ylabel", metavar="TEXT",
                   help="Y axis title (default: column name with units)")
    p.add_argument("--no-legend", action="store_true", help="Do not draw a legend")
    p.add_argument("--base-size", type=float, default=11, metavar="PT",
                   help="Base font size in points (default: %(default)s)")
    p.add_argument("--point-size", type=float, default=25,
                   help="Marker size (default: %(default)s)")
    p.add_argument("--width", type=float, default=7, help="Figure width in inches (default: %(default)s)")
    p.add_argument("--height", type=float, default=5, help="Figure height in inches (default: %(default)s)")
    return p, p.parse_args()


def read_bin_stats(path, parser):
    """Load a CheckM bin_stats*.tsv (name<TAB>python-dict per line) as a table."""
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            name, sep, stats = line.partition("\t")
            try:
                stats = ast.literal_eval(stats) if sep else None
            except (ValueError, SyntaxError):
                stats = None
            if not isinstance(stats, dict):
                parser.error(f"{path} line {lineno}: expected '<bin name><TAB>{{...}}' "
                             "(an official CheckM bin_stats file)")
            row = {NAME_COL: name}
            for key, value in stats.items():
                # GCN0..GCN5+ hold lists of marker IDs, which cannot be plotted
                if not isinstance(value, (list, tuple, dict, set)):
                    row[RENAME.get(key, key)] = value
            rows.append(row)
    if not rows:
        parser.error(f"input file is empty: {path}")
    df = pd.DataFrame(rows)
    if df[NAME_COL].duplicated().any():
        dup = df.loc[df[NAME_COL].duplicated(), NAME_COL].iloc[0]
        print(f"Warning: duplicate genome names (e.g. {dup})", file=sys.stderr)
    return df


def resolve_column(df, spec, parser, flag):
    """Return the real column name for a user-supplied name or 1-based index."""
    if spec in df.columns and not spec.isdigit():
        return spec
    if spec.isdigit():
        idx = int(spec)
        if 1 <= idx <= len(df.columns):
            return df.columns[idx - 1]
        parser.error(f"{flag}: column number {idx} out of range (1-{len(df.columns)})")
    lookup = {c.lower().strip(): c for c in df.columns}
    key = spec.lower().strip()
    if key in lookup:
        return lookup[key]
    hint = ""
    if key in ("completeness", "contamination") or key in {v.lower() for v in RENAME.values()}:
        hint = (f" ('{spec}' is only in bin_stats_ext.tsv; bin_stats.analyze.tsv and "
                "bin_stats.tree.tsv hold assembly stats only)")
    parser.error(f"{flag}: unknown column '{spec}'{hint}. Available: "
                 + ", ".join(f"'{c}'" for c in df.columns))


def to_numeric(df, col, parser, flag):
    values = pd.to_numeric(df[col], errors="coerce")
    if values.notna().sum() == 0:
        parser.error(f"{flag}: column '{col}' is not numeric and cannot be plotted")
    return values


def padded_limits(lo, hi, frac, log, parser, flag):
    """Widen [lo, hi] by frac of the span on each side (in log space for log axes)."""
    if lo >= hi:
        parser.error(f"{flag}: MIN must be smaller than MAX (got {lo:g} {hi:g})")
    if log:
        if lo <= 0:
            parser.error(f"{flag}: limits must be > 0 on a log-scaled axis")
        llo, lhi = math.log10(lo), math.log10(hi)
        pad = (lhi - llo) * frac
        return 10 ** (llo - pad), 10 ** (lhi + pad)
    pad = (hi - lo) * frac
    return lo - pad, hi + pad


def read_name_list(path, parser):
    if not os.path.isfile(path):
        parser.error(f"--exclude: file not found: {path}")
    names = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                names.add(line.split()[0])
    return names


def failing_mask(values, thresholds, rule):
    """True where a genome fails 'value >= / <= strictest threshold' (NaN fails)."""
    if rule == "ge":
        cut = max(thresholds)
        return ~(values >= cut), cut
    cut = min(thresholds)
    return ~(values <= cut), cut


def safe_name(text):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9.-]+", "_", text)).strip("_")


def apply_theme(base):
    """Approximate ggplot2::theme_grey(base_size=base), the ggplot default."""
    tick = base * 0.8
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": FONTS,
        "font.size": base,
        "mathtext.default": "regular",
        "axes.facecolor": GREY92,
        "axes.edgecolor": "none",
        "axes.linewidth": 0,
        "axes.spines.left": False,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.spines.bottom": False,
        "axes.labelsize": base,
        "axes.labelcolor": "black",
        "axes.labelpad": 5,
        "axes.titlesize": base * 1.2,
        "axes.titlelocation": "left",
        "axes.titlepad": 8,
        "axes.grid": True,
        "axes.grid.which": "both",
        "axes.axisbelow": True,
        "grid.color": "white",
        "grid.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.color": GREY20,
        "ytick.color": GREY20,
        "xtick.labelcolor": GREY30,
        "ytick.labelcolor": GREY30,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.size": 0,
        "ytick.minor.size": 0,
        "legend.frameon": False,
        "legend.fontsize": tick,
        "legend.title_fontsize": base * 0.9,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def axis_title(col, values, override):
    """Return (axis title, display divisor) for a column."""
    if override:
        return override, 1
    label, unit = PRETTY.get(col, (col, None))
    scale = 1
    if unit == "bp":
        top = values.abs().max()
        if top >= 1e6:
            unit, scale = "Mbp", 1e6
        elif top >= 1e4:
            unit, scale = "kbp", 1e3
    return (f"{label} ({unit})" if unit else label), scale


class TickFormatter(Formatter):
    """ggplot-like tick labels: same decimals on every tick, display scaling
    (e.g. bp -> Mbp), and thousands separators once values reach 10,000."""

    def __init__(self, scale=1, log=False):
        self.scale = scale
        self.log = log
        self.decimals = 0
        self.commas = False

    def set_locs(self, locs):
        super().set_locs(locs)
        vals = [abs(v) / self.scale for v in locs]
        self.commas = bool(vals) and max(vals) >= 1e4
        self.decimals = max((len(f"{v:.6f}".rstrip("0").split(".")[1]) for v in vals),
                            default=0)

    def __call__(self, value, pos=None):
        v = value / self.scale
        if self.log:
            # log ticks span decades, so give each its own minimal precision
            return f"{v:,.0f}" if abs(v) >= 1e4 else f"{v:g}"
        return f"{v:{',' if self.commas else ''}.{self.decimals}f}"


def style_axis(axis, scale, log, lo, hi):
    axis.set_major_formatter(TickFormatter(scale, log))
    if not log:
        axis.set_minor_locator(AutoMinorLocator(2))
        axis.set_minor_formatter(NullFormatter())
    elif hi / lo < 30:
        # less than ~1.5 decades: label the minor ticks too, otherwise none show
        axis.set_minor_formatter(TickFormatter(scale, log))
    else:
        axis.set_minor_formatter(NullFormatter())


def place_legend_below(fig, ax, handles, labels, gap_pt=6):
    """Centre the legend under the x-axis title, as wide a row as fits the panel."""
    renderer = fig.canvas.get_renderer()
    bottom_px = ax.get_tightbbox(renderer).y0 - gap_pt * fig.dpi / 72
    y = ax.transAxes.inverted().transform((0, bottom_px))[1]
    panel_w = ax.get_window_extent(renderer).width
    for ncol in range(len(handles), 0, -1):
        leg = ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y),
                        ncol=ncol, borderaxespad=0, handlelength=1.2,
                        handletextpad=0.5, columnspacing=1.5)
        if ncol == 1 or leg.get_window_extent(renderer).width <= panel_w:
            return leg
        leg.remove()


def threshold_failures(x, y, args, xcol, ycol):
    """Mask of genomes failing --xline/--yline, plus file-name tags and rule text."""
    fail = pd.Series(False, index=x.index)
    tags, rules = [], []
    if args.xline:
        fx, cut = failing_mask(x, args.xline, args.x_pass)
        fail |= fx
        tags.append(f"{safe_name(xcol)}_{args.x_pass}{cut:g}")
        rules.append(f"{xcol} {'>=' if args.x_pass == 'ge' else '<='} {cut:g}")
    if args.yline:
        fy, cut = failing_mask(y, args.yline, args.y_pass)
        fail |= fy
        tags.append(f"{safe_name(ycol)}_{args.y_pass}{cut:g}")
        rules.append(f"{ycol} {'>=' if args.y_pass == 'ge' else '<='} {cut:g}")
    return fail, tags, rules


def repel_labels(fig, ax, lx, ly, names, all_x, all_y, point_size, fontsize):
    """ggrepel-like labelling: put each name at the nearest spot (searching
    outward around its point) where the label box overlaps no data point, no
    other label and stays inside the panel; connect it with a leader line."""
    renderer = fig.canvas.get_renderer()
    panel = ax.get_window_extent(renderer)
    px_per_pt = fig.dpi / 72
    radius = math.sqrt(point_size) / 2 * px_per_pt   # scatter s is marker area in pt^2
    gap = 2 * px_per_pt                              # clearance around points/labels
    box_pad = 0.25                                   # label box padding (x fontsize)
    to_disp = ax.transData.transform
    to_data = ax.transData.inverted().transform
    pts = to_disp(np.column_stack([all_x, all_y]))
    targets = to_disp(np.column_stack([lx, ly]))

    inside = ((targets[:, 0] >= panel.x0) & (targets[:, 0] <= panel.x1)
              & (targets[:, 1] >= panel.y0) & (targets[:, 1] <= panel.y1))
    items = []
    for (tx, ty), name, ok in zip(targets, names, inside):
        if not ok:
            continue  # point is outside the zoomed axis range
        txt = ax.text(0, 0, name, fontsize=fontsize, color="black", ha="center",
                      va="center", zorder=5,
                      bbox=dict(boxstyle=f"round,pad={box_pad}", fc="white",
                                ec=GREY30, lw=0.4, alpha=0.9))
        ext = txt.get_window_extent(renderer)
        pad_px = 2 * box_pad * fontsize * px_per_pt
        items.append([txt, tx, ty, ext.width + pad_px, ext.height + pad_px])

    # most crowded points pick first, while there is still free space nearby
    def crowding(item):
        d = np.hypot(pts[:, 0] - item[1], pts[:, 1] - item[2])
        return -int((d < 40 * px_per_pt).sum())
    items.sort(key=crowding)

    # preferred directions first (upper-right, right, up, ...), then the rest
    angles = [45, 0, 90, 135, 180, -45, -90, -135, 22.5, 67.5, 112.5, 157.5,
              -22.5, -67.5, -112.5, -157.5]
    dirs = [(math.cos(math.radians(a)), math.sin(math.radians(a))) for a in angles]
    step = 4 * px_per_pt
    placed = []
    crowded = []
    def overlaps(cx, cy, w, h):
        x0, x1, y0, y1 = cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2
        if x0 < panel.x0 + gap or x1 > panel.x1 - gap or y0 < panel.y0 + gap or y1 > panel.y1 - gap:
            return 1e9
        # distance from every point centre to the box
        dx = np.maximum(np.maximum(x0 - pts[:, 0], 0), pts[:, 0] - x1)
        dy = np.maximum(np.maximum(y0 - pts[:, 1], 0), pts[:, 1] - y1)
        hits = int((np.hypot(dx, dy) < radius + gap).sum())
        for bx0, bx1, by0, by1 in placed:
            if x0 < bx1 + gap and x1 > bx0 - gap and y0 < by1 + gap and y1 > by0 - gap:
                hits += 10
        return hits

    for txt, tx, ty, w, h in items:
        best = None
        for k in range(60):
            dist = radius + gap + k * step
            for cos_a, sin_a in dirs:
                cx = tx + cos_a * (dist + w / 2)
                cy = ty + sin_a * (dist + h / 2)
                score = overlaps(cx, cy, w, h)
                if best is None or score < best[0]:
                    best = (score, cx, cy)
                if score == 0:
                    break
            if best[0] == 0:
                break
        score, cx, cy = best
        if score:
            crowded.append(txt.get_text())
        placed.append((cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2))
        txt.set_position(to_data((cx, cy)))

        # leader line from the marker edge to the nearest point of the label box
        ex = min(max(tx, cx - w / 2), cx + w / 2)
        ey = min(max(ty, cy - h / 2), cy + h / 2)
        length = math.hypot(ex - tx, ey - ty)
        if length > radius + 3 * px_per_pt:
            sx = tx + (ex - tx) / length * radius
            sy = ty + (ey - ty) / length * radius
            (x0d, y0d), (x1d, y1d) = to_data((sx, sy)), to_data((ex, ey))
            ax.plot([x0d, x1d], [y0d, y1d], color=GREY20, linewidth=0.5,
                    zorder=4, solid_capstyle="round")
    if crowded:
        print(f"Warning: no overlap-free spot for {len(crowded)} label(s) "
              f"(e.g. {crowded[0]}); try a larger --width/--height, a smaller "
              "--base-size or --point-size, or zoom with --xlim/--ylim", file=sys.stderr)
    return len(items)


def strictest_first(values, rule):
    return sorted(values, reverse=(rule == "ge"))


def shade_failing(ax, xs, ys, x_pass, y_pass, color, alpha, labels):
    """Tint everything outside the passing rectangle (one tier per threshold).

    Returns legend proxy handles whose opacity matches each stacked tier."""
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    xs = strictest_first(xs or [], x_pass)
    ys = strictest_first(ys or [], y_pass)
    outer = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    handles = []
    for k in range(max(len(xs), len(ys))):
        xc = xs[min(k, len(xs) - 1)] if xs else None
        yc = ys[min(k, len(ys) - 1)] if ys else None
        px = (x0, x1) if xc is None else ((xc, x1) if x_pass == "ge" else (x0, xc))
        py = (y0, y1) if yc is None else ((yc, y1) if y_pass == "ge" else (y0, yc))
        px = (max(px[0], x0), min(px[1], x1))
        py = (max(py[0], y0), min(py[1], y1))
        verts = list(outer)
        codes = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
        if px[0] < px[1] and py[0] < py[1]:
            # opposite winding cuts the passing rectangle out as a hole
            verts += [(px[0], py[0]), (px[0], py[1]), (px[1], py[1]),
                      (px[1], py[0]), (px[0], py[0])]
            codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, alpha=alpha,
                               edgecolor="none", zorder=0.8))
        a = 1 - (1 - alpha) ** (k + 1)
        # pre-blend over the grey panel so the legend swatch matches the plot
        swatch = tuple(a * c + (1 - a) * g for c, g in zip(to_rgb(color), to_rgb(GREY92)))
        handles.append(Patch(facecolor=swatch, edgecolor="none", label=labels[k]))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    return handles


def fail_rule_texts(xcol, ycol, xs, ys, x_pass, y_pass):
    """One 'X < a or Y > b' description per shading tier, strictest first."""
    xs = strictest_first(xs or [], x_pass)
    ys = strictest_first(ys or [], y_pass)
    texts = []
    for k in range(max(len(xs), len(ys))):
        parts = []
        if xs:
            parts.append(f"{xcol} {'<' if x_pass == 'ge' else '>'} {xs[min(k, len(xs) - 1)]:g}")
        if ys:
            parts.append(f"{ycol} {'<' if y_pass == 'ge' else '>'} {ys[min(k, len(ys) - 1)]:g}")
        texts.append(" or ".join(parts))
    return texts


def main():
    parser, args = parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"input file not found: {args.input}")
    df = read_bin_stats(args.input, parser)
    name_col = NAME_COL

    if args.list:
        print(f"Columns in {args.input} ({len(df)} genomes):")
        for i, c in enumerate(df.columns, 1):
            numeric = c != name_col and pd.to_numeric(df[c], errors="coerce").notna().any()
            note = "" if numeric or c == name_col else "   (text - not plottable)"
            print(f"  {i:>2}  {c}{note}")
        return

    if not args.xcol or not args.ycol:
        parser.error("both -x and -y are required (use --list to see column choices)")
    if args.write_failed and not (args.xline or args.yline):
        parser.error("--write-failed needs at least one --xline or --yline threshold")
    if args.label_failed and not (args.xline or args.yline or args.highlight_excluded):
        parser.error("--label-failed needs --xline/--yline thresholds and/or "
                     "--exclude FILE --highlight-excluded")
    if args.highlight_excluded and not args.exclude:
        parser.error("--highlight-excluded requires at least one --exclude FILE")

    highlighted = pd.Series(False, index=df.index)
    if args.exclude:
        excluded = set()
        for path in args.exclude:
            excluded |= read_name_list(path, parser)
        names = df[name_col].astype(str)
        hit = names.isin(excluded)
        missing = excluded - set(names)
        action = "Retained for highlighting" if args.highlight_excluded else "Excluded"
        print(f"{action} {int(hit.sum())} genome(s) listed in: {', '.join(args.exclude)}")
        if missing:
            print(f"Warning: {len(missing)} name(s) in the exclude list were not found "
                  f"in {args.input} (e.g. {sorted(missing)[0]})", file=sys.stderr)
        if args.highlight_excluded:
            highlighted = hit
        else:
            df = df[~hit].reset_index(drop=True)
            highlighted = pd.Series(False, index=df.index)
        if df.empty:
            parser.error("no genomes left to plot after --exclude")

    xcol = resolve_column(df, args.xcol, parser, "-x")
    ycol = resolve_column(df, args.ycol, parser, "-y")
    x = to_numeric(df, xcol, parser, "-x")
    y = to_numeric(df, ycol, parser, "-y")

    fail = None
    if args.xline or args.yline:
        fail, tags, rules = threshold_failures(x, y, args, xcol, ycol)

    if args.write_failed:
        fail_out = (args.write_failed if args.write_failed != AUTO
                    else f"CheckM_failed_{'_'.join(tags)}.txt")
        failed = df.loc[fail, name_col].astype(str)
        with open(fail_out, "w") as fh:
            fh.writelines(f"{n}\n" for n in failed)
        print(f"Pass criteria: {' AND '.join(rules)}")
        print(f"{len(failed)} of {len(df)} genome(s) FAILED -> written to {fail_out}")

    c = None
    ccol = None
    if args.color_by:
        ccol = resolve_column(df, args.color_by, parser, "--color-by")
        c = to_numeric(df, ccol, parser, "--color-by")

    mask = x.notna() & y.notna()
    if c is not None:
        mask &= c.notna()
    dropped = int((~mask).sum())
    if dropped:
        print(f"Warning: skipping {dropped} row(s) with missing/non-numeric values", file=sys.stderr)

    apply_theme(args.base_size)
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    regular = mask & ~highlighted
    marked = mask & highlighted
    sc = ax.scatter(x[regular], y[regular], s=args.point_size,
                    c=c[regular] if c is not None else POINT_COLOR,
                    cmap="viridis" if c is not None else None,
                    edgecolors=GREY20, linewidths=0.35, alpha=0.85, zorder=2,
                    label="QA Pass" if args.highlight_excluded and regular.any() else None)
    cbar = None
    if c is not None and regular.any():
        cbar = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03, aspect=25)
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(length=2, width=0.5)
        cbar.set_label(axis_title(ccol, c[mask], None)[0])
        cbar.ax.yaxis.set_major_formatter(TickFormatter())
    if args.highlight_excluded:
        if marked.any():
            ax.scatter(x[marked], y[marked], s=args.point_size * 1.1,
                       color=HIGHLIGHT_COLOR, marker="D", edgecolors=GREY20,
                       linewidths=0.5, zorder=3, label="QA Pass")
        print(f"Plotted {int(marked.sum())} highlighted genome(s) and "
              f"{int(regular.sum())} other genome(s)")

    if args.label:
        for name, xv, yv in zip(df.loc[mask, name_col], x[mask], y[mask]):
            ax.annotate(str(name), (xv, yv), fontsize=args.base_size * 0.45,
                        color=GREY20, xytext=(2, 2), textcoords="offset points",
                        zorder=4)

    if args.logx:
        ax.set_xscale("log")
    if args.logy:
        ax.set_yscale("log")

    # make sure the thresholds (and so the shaded region) are inside the view
    if (args.xline or args.yline) and mask.any():
        xm, ym = float(x[mask].median()), float(y[mask].median())
        pts = [(v, ym) for v in args.xline or [] if not (args.logx and v <= 0)]
        pts += [(xm, v) for v in args.yline or [] if not (args.logy and v <= 0)]
        if pts:
            ax.update_datalim(pts)
            ax.autoscale_view()

    if args.xlim:
        ax.set_xlim(*padded_limits(*args.xlim, args.pad, args.logx, parser, "--xlim"))
    if args.ylim:
        ax.set_ylim(*padded_limits(*args.ylim, args.pad, args.logy, parser, "--ylim"))

    shade_handles = []
    if args.xline or args.yline:
        if not 0 <= args.shade_alpha <= 1:
            parser.error("--shade-alpha must be between 0 and 1")
        tier_labels = fail_rule_texts(PRETTY.get(xcol, (xcol,))[0], PRETTY.get(ycol, (ycol,))[0],
                                      args.xline, args.yline, args.x_pass, args.y_pass)
        if args.shade_label:
            tier_labels[0] = args.shade_label
        shade_handles = shade_failing(ax, args.xline, args.yline, args.x_pass, args.y_pass,
                                      args.shade_color, args.shade_alpha, tier_labels)

    xlabel, xscale = axis_title(xcol, x[mask], args.xlabel)
    ylabel, yscale = axis_title(ycol, y[mask], args.ylabel)
    style_axis(ax.xaxis, xscale, args.logx, *ax.get_xlim())
    style_axis(ax.yaxis, yscale, args.logy, *ax.get_ylim())
    ax.grid(which="major", color="white", linewidth=0.8, linestyle="-")
    ax.grid(which="minor", color="white", linewidth=0.5, linestyle=(0, (3, 2)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if not args.no_title:
        xname, yname = PRETTY.get(xcol, (xcol,))[0], PRETTY.get(ycol, (ycol,))[0]
        ax.set_title(args.title or f"{yname} vs {xname} (n={int(mask.sum())})")

    if args.label_failed:
        lab = mask & highlighted
        if fail is not None:
            lab |= mask & fail
        n = repel_labels(fig, ax, x[lab].to_numpy(), y[lab].to_numpy(),
                         df.loc[lab, name_col].astype(str).tolist(),
                         x[mask].to_numpy(), y[mask].to_numpy(),
                         args.point_size, args.base_size * 0.65)
        what = []
        if fail is not None:
            what.append("in the failing (red) region")
        if args.highlight_excluded:
            what.append("highlighted from --exclude")
        print(f"Labelled {n} genome(s) {' or '.join(what)}")

    handles, labels = ax.get_legend_handles_labels()
    handles += shade_handles
    labels += [h.get_label() for h in shade_handles]
    if handles and not args.no_legend:
        place_legend_below(fig, ax, handles, labels)

    out = args.output or f"CheckM_{safe_name(ycol)}_vs_{safe_name(xcol)}.png"
    if not out.lower().endswith(".png"):
        out += ".png"
    fig.savefig(out, dpi=300, format="png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Saved plot: {out}")


if __name__ == "__main__":
    main()
