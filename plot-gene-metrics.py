#!/usr/bin/env python3

# Plot the distribution of 7 PhyKIT information-content metrics across all
# single-copy BUSCO genes as violin + box plots, one panel per metric.
# Input is the tab-delimited output from compute_info_content.py.
#
# Reference: https://jlsteenwyk.com/tutorials/ub_sensitivity_analysis.html#Section2
#
# Dependencies:
#   - matplotlib
#   - pandas

import os
import sys
import argparse

from time import gmtime, strftime

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


METRIC_ORDER = [
    "aln_len",
    "abs",
    "rcv",
    "lbs",
    "treeness",
    "saturation",
    "treeness_over_rcv",
]

METRIC_LABELS = {
    "aln_len":           "Alignment length",
    "abs":               "Avg. bipartition support",
    "rcv":               "Relative composition\nvariability (RCV)",
    "lbs":               "Median long-branch\nscore (LBS)",
    "treeness":          "Treeness",
    "saturation":        "Saturation",
    "treeness_over_rcv": "Treeness / RCV",
}

# Metrics where lower values indicate better phylogenetic signal
LOWER_IS_BETTER = {"rcv", "lbs", "saturation"}


def main():
    parser = argparse.ArgumentParser(
        description="Plot PhyKIT information-content metrics from compute_info_content.py output"
    )
    parser.add_argument(
        "--input", type=str, required=True,
        help="Tab-delimited input file produced by compute_info_content.py (gene, metric, value)"
    )
    parser.add_argument(
        "--output", type=str, default="info_content_genes.pdf",
        help="Output figure file (default: info_content_genes.pdf). Format inferred from extension (pdf, png, svg)."
    )
    parser.add_argument(
        "--dpi", type=int, default=150,
        help="Resolution in DPI for raster formats (default: 150)"
    )
    args = parser.parse_args()

    input_file = os.path.abspath(args.input)
    output_file = os.path.abspath(args.output)
    dpi = args.dpi

    if not os.path.isfile(input_file):
        print("Error! " + input_file + " does not exist!")
        sys.exit(1)

    # Load data
    print_message("Reading " + input_file)
    df = pd.read_csv(input_file, sep="\t", header=None, names=["gene", "metric", "value"])
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value"])

    metrics_found = df["metric"].unique().tolist()
    metrics_to_plot = [m for m in METRIC_ORDER if m in metrics_found]
    missing = [m for m in METRIC_ORDER if m not in metrics_found]

    if missing:
        print_message("Warning: metrics not found in input and will be skipped: " + ", ".join(missing))
    if not metrics_to_plot:
        print_message("Error! No known metrics found in input file. Exiting.")
        sys.exit(1)

    n_genes = df["gene"].nunique()
    print_message(str(n_genes) + " genes and " + str(len(metrics_to_plot)) + " metrics loaded")

    # Layout: up to 4 columns
    n_cols = 4
    n_rows = (len(metrics_to_plot) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 4.5 * n_rows))
    axes = axes.flatten()

    for i, metric in enumerate(metrics_to_plot):
        ax = axes[i]
        values = df[df["metric"] == metric]["value"].dropna().tolist()

        # Violin
        parts = ax.violinplot(values, positions=[0], showmedians=False, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor("#4C72B0")
            body.set_alpha(0.6)
            body.set_edgecolor("#2a4a80")
            body.set_linewidth(0.8)

        # Box plot overlaid
        bp = ax.boxplot(
            values,
            positions=[0],
            widths=0.12,
            patch_artist=True,
            boxprops=dict(facecolor="white", linewidth=1.2),
            medianprops=dict(color="#C44E52", linewidth=2),
            whiskerprops=dict(linewidth=1.2),
            capprops=dict(linewidth=1.2),
            flierprops=dict(marker="o", markerfacecolor="#4C72B0", markersize=3, alpha=0.5, linewidth=0),
        )

        # Annotate with n, median, mean
        median = pd.Series(values).median()
        mean = pd.Series(values).mean()
        ax.text(
            0.97, 0.97,
            f"n = {len(values)}\nmedian = {median:.3f}\nmean = {mean:.3f}",
            transform=ax.transAxes,
            fontsize=8,
            va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="lightgrey", alpha=0.8),
        )

        # Signal direction hint
        direction = "lower = better" if metric in LOWER_IS_BETTER else "higher = better"
        ax.text(
            0.03, 0.03,
            direction,
            transform=ax.transAxes,
            fontsize=7,
            va="bottom", ha="left",
            color="grey",
            style="italic",
        )

        ax.set_title(METRIC_LABELS.get(metric, metric), fontsize=11, fontweight="bold", pad=8)
        ax.set_xticks([])
        ax.set_ylabel("Value", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

    # Hide unused panels
    for j in range(len(metrics_to_plot), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "PhyKIT information-content metrics across single-copy BUSCO genes",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()

    print_message("Saving figure to " + output_file)
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    print_message("Done!")


def print_message(*message):
    print(strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message)))


if __name__ == "__main__":
    main()
