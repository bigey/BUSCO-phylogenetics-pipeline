#!/usr/bin/env python3

# Phylogenomic sensitivity analysis: for each PhyKIT information-content metric,
# select the top-scoring genes, concatenate their alignments using PhyKIT
# create_concat, and infer a phylogenetic tree with IQ-TREE.
# Optionally computes Robinson-Foulds distances against a reference tree.
#
# Input is the tab-delimited output from compute-gene-metrics.py,
#    which should have three columns: gene, metric, value.
# Reference: https://jlsteenwyk.com/tutorials/ub_sensitivity_analysis.html#Section2
#
# Dependencies:
#   - PhyKIT
#   - IQ-TREE

import os
import sys
import shutil
import argparse
import logging
import colorlog
import subprocess
import multiprocessing as mp


# Metrics where lower values indicate better phylogenetic signal
LOWER_IS_BETTER = {"rcv", "lbs", "saturation"}

METRIC_ORDER = [
    "aln_len",
    "abs",
    "rcv",
    "lbs",
    "treeness",
    "saturation",
    "treeness_over_rcv",
]

external_dependencies = ["phykit", "iqtree3"]


def check_dependency(name):
    if shutil.which(name) is None:
        logger.critical(f"Program {name} is required but not in the PATH!")
        logger.critical(f"Please install {name}.")
        sys.exit(1)


def main(args):
    input_file = os.path.abspath(args.input)
    aln_dir = os.path.abspath(args.trimmed_alignments)
    output_dir = os.path.abspath(args.output_dir)
    top_fraction = args.top_fraction
    model = args.model
    threads = args.threads
    ref_tree = os.path.abspath(args.ref_tree) if args.ref_tree else None

    # Validate inputs
    if not os.path.isfile(input_file):
        logger.critical(f"{input_file} does not exist!")
        sys.exit(1)

    if not os.path.isdir(aln_dir):
        logger.critical(f"{aln_dir} is not a directory!")
        sys.exit(1)

    if os.path.isdir(output_dir):
        logger.critical(f"{output_dir} already exists!")
        sys.exit(1)

    if ref_tree and not os.path.isfile(ref_tree):
        logger.critical(f"Reference tree {ref_tree} does not exist!")
        sys.exit(1)

    if not (0 < top_fraction <= 1):
        logger.critical("--top-fraction must be between 0 (exclusive) and 1 (inclusive)!")
        sys.exit(1)

    os.mkdir(output_dir)

    logger.info(f"Reading {input_file}")
    genes_by_metric = {}

    with open(input_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 3:
                continue
            gene, metric, value = parts
            if metric not in METRIC_ORDER:
                continue
            try:
                value = float(value)
            except ValueError:
                continue
            if metric not in genes_by_metric:
                genes_by_metric[metric] = []
            genes_by_metric[metric].append((gene, value))

    if not genes_by_metric:
        logger.critical(f"No valid metric data found in {input_file}. Exiting.")
        sys.exit(1)

    metrics_found = [m for m in METRIC_ORDER if m in genes_by_metric]
    logger.info(f"{len(metrics_found)} metrics loaded: {', '.join(metrics_found)}")

    # Split threads between parallel workers and IQ-TREE so total CPU usage
    # stays close to --threads (e.g. 64 threads // 7 metrics = 9 threads per IQ-TREE)
    n_metrics = len(metrics_found)
    iqtree_threads = max(1, threads // n_metrics)
    pool_size = min(threads, n_metrics)

    jobs = []
    n_top = 1
    for metric in metrics_found:
        gene_values = genes_by_metric[metric]
        n_top = max(1, int(len(gene_values) * top_fraction))
        jobs.append((metric, gene_values, aln_dir, output_dir, n_top, model, ref_tree, iqtree_threads))

    logger.info(f"Running sensitivity analysis of the best scoring genes (top {top_fraction:.2f} fraction)")
    logger.info(f"Running {pool_size} parallel workers, {iqtree_threads} IQ-TREE thread(s) each")

    pool = mp.Pool(processes=pool_size)
    results = pool.map(run_metric, jobs)
    pool.close()
    pool.join()

    # Write RF distance summary if reference tree was provided
    if ref_tree:
        rf_file = os.path.join(output_dir, "rf-distances.tsv")
        logger.info(f"Writing Robinson-Foulds distances to {rf_file}")
        with open(rf_file, "w") as fout:
            fout.write("metric\tt_genes\tgenes\trf_distance\trf_normalized\n")
            for row in results:
                fout.write(f"{row['metric']}\t{row['t_genes']}\t{row['genes']}\t{row['rf']}\t{row['norm_rf']}\n")

    logger.info(f"All results are in {output_dir}")
    logger.info(f"Sensitivity analysis complete!")
    logger.info(f"All done. Exit")


def run_metric(args):
    metric, gene_values, aln_dir, output_dir, n_top, model, ref_tree, iqtree_threads = args

    # Sort: ascending for metrics where lower = better, descending otherwise
    t_genes = len(gene_values)
    descending = metric in LOWER_IS_BETTER
    sorted_genes = sorted(gene_values, key=lambda x: x[1], reverse=not descending)
    top_genes = [gene for gene, _ in sorted_genes[:n_top]]

    logger.info(f"{metric}: selected {n_top} / {len(gene_values)} genes")

    # Write gene list file (one alignment path per line)
    gene_list_file = os.path.join(output_dir, f"{metric}.gene_list.txt")
    with open(gene_list_file, "w") as fout:
        for gene in top_genes:
            fout.write(os.path.join(aln_dir, f"{gene}.trimmed.aln.fasta") + "\n")

    # Concatenate alignments with PhyKIT
    prefix = os.path.join(output_dir, metric)
    ret = subprocess.run(
        ["phykit", "create_concat", "-a", gene_list_file, "-p", prefix],
        capture_output=True, text=True
    )
    concat_fasta = f"{prefix}.fa"
    if ret.returncode != 0 or not os.path.isfile(concat_fasta):
        logger.warning(f"phykit create_concat failed for metric {metric}")
        return {"metric": metric, "t_genes": t_genes, "genes": n_top, "rf": "NA", "norm_rf": "NA"}

    logger.debug(f"{metric}: alignment concatenated to {concat_fasta}")

    # Infer phylogenetic tree with IQ-TREE
    subprocess.run(
        ["iqtree3", "--quiet", "-s", concat_fasta, "-m", model, "--fast", "-T", str(iqtree_threads)],
        capture_output=True
    )

    treefile = f"{concat_fasta}.treefile"
    if not os.path.isfile(treefile):
        logger.warning(f"IQ-TREE did not produce a treefile for metric {metric}")
        return {"metric": metric, "t_genes": t_genes, "genes": n_top, "rf": "NA", "norm_rf": "NA"}

    logger.debug(f"{metric}: tree inferred at {treefile}")

    plain_rf = "NA"
    norma_rf = "NA"

    # Compute Robinson-Foulds distance against reference tree
    if ref_tree:
        out = run_command(["phykit", "robinson_foulds_distance", ref_tree, treefile])
        if out:
            plain_rf = out.strip().split()[0]
            norma_rf = out.strip().split()[1]
            logger.debug(f"{metric}: RF distance to reference = {plain_rf} / {norma_rf}")
        else:
            logger.warning(f"Could not compute RF distance for metric {metric}")

    return {"metric": metric, "t_genes": t_genes, "genes": n_top, "rf": plain_rf, "norm_rf": norma_rf}


def run_command(cmd):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return None
        return result.stdout
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Phylogenomic sensitivity analysis using PhyKIT information-content metrics"
    )
    parser.add_argument(
        "--input", type=str, required=True,
        help="Tab-delimited (TSV) file produced by compute-gene-metrics.py (gene, metric, value)"
    )
    parser.add_argument(
        "--trimmed-alignments", type=str, required=True,
        help="Directory containing trimmed alignment files (*.trimmed.aln.fasta)"
    )
    parser.add_argument(
        "--output-dir", type=str, required=True,
        help="Output directory to store results (must not already exist)"
    )
    parser.add_argument(
        "--top-fraction", type=float, default=0.75,
        help="Fraction of top-scoring genes to retain per metric (default: 0.75 = top 75%%)"
    )
    parser.add_argument(
        "--model", type=str, default="LG+R4+F",
        help="IQ-TREE substitution model (default: LG+R4+F)"
    )
    parser.add_argument(
        "--threads", type=int, default=8,
        help="Number of parallel threads for running metrics concurrently (default: 8)"
    )
    parser.add_argument(
        "--ref-tree", type=str, default=None,
        help="Optional reference tree to compute Robinson-Foulds distance against each sensitivity tree"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Turn on verbose mode."
    )
    args = parser.parse_args()

    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    log_colors = {
        "DEBUG": "cyan",
        "INFO": "green",
        "WARNING": "yellow",
        "ERROR": "white,bg_red",
        "CRITICAL": "red",
    }

    formatter = colorlog.ColoredFormatter(
        fmt="%(asctime)s:%(log_color)s%(levelname)s%(reset)s:%(message)s",
        log_colors=log_colors,
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    handler = colorlog.StreamHandler()
    handler.setFormatter(fmt=formatter)

    logger = logging.getLogger()
    logger.addHandler(handler)
    logger.setLevel(log_level)

    list(map(check_dependency, external_dependencies))

    main(args)
