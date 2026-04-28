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
import argparse
import subprocess
import multiprocessing as mp

from time import gmtime, strftime


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


def main():
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
        "--top-fraction", type=float, default=0.9,
        help="Fraction of top-scoring genes to retain per metric (default: 0.9 = top 90%%)"
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
    args = parser.parse_args()

    input_file = os.path.abspath(args.input)
    aln_dir = os.path.abspath(args.trimmed_alignments)
    output_dir = os.path.abspath(args.output)
    top_fraction = args.top_fraction
    model = args.model
    threads = args.threads
    ref_tree = os.path.abspath(args.ref_tree) if args.ref_tree else None

    # Validate inputs
    if not os.path.isfile(input_file):
        print("Error! " + input_file + " does not exist!")
        sys.exit(1)

    if not os.path.isdir(aln_dir):
        print("Error! " + aln_dir + " is not a directory!")
        sys.exit(1)

    if os.path.isdir(output_dir):
        print("Error! " + output_dir + " already exists!")
        sys.exit(1)

    if ref_tree and not os.path.isfile(ref_tree):
        print("Error! Reference tree " + ref_tree + " does not exist!")
        sys.exit(1)

    if not (0 < top_fraction <= 1):
        print("Error! --top-percent must be between 0 (exclusive) and 1 (inclusive)!")
        sys.exit(1)

    os.mkdir(output_dir)

    # Load metrics TSV
    print_message("Reading " + input_file)
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
        print_message("Error! No valid metric data found in " + input_file + ". Exiting.")
        sys.exit(1)

    metrics_found = [m for m in METRIC_ORDER if m in genes_by_metric]
    print_message(str(len(metrics_found)) + " metrics loaded: " + ", ".join(metrics_found))

    # Split threads between parallel workers and IQ-TREE so total CPU usage
    # stays close to --threads (e.g. 64 threads / 7 metrics = 9 threads per IQ-TREE)
    n_metrics = len(metrics_found)
    iqtree_threads = max(1, threads // n_metrics)
    pool_size = min(threads, n_metrics)

    jobs = []
    for metric in metrics_found:
        gene_values = genes_by_metric[metric]
        n_top = max(1, int(len(gene_values) * top_fraction))
        jobs.append((metric, gene_values, aln_dir, output_dir, n_top, model, ref_tree, iqtree_threads))

    pct_label = str(int(top_fraction * 100)) + "%"
    print_message(
        "Running sensitivity analysis for " + str(n_metrics) + " metrics "
        "(top " + pct_label + " genes per metric) using " + str(pool_size) + " parallel workers, "
        + str(iqtree_threads) + " IQ-TREE thread(s) each"
    )

    pool = mp.Pool(processes=pool_size)
    results = pool.map(run_metric, jobs)
    pool.close()
    pool.join()

    # Print per-metric summaries
    for log_lines, _ in results:
        for line in log_lines:
            print(line)

    # Write RF distance summary if reference tree was provided
    if ref_tree:
        rf_file = os.path.join(output_dir, "rf_distances.tsv")
        print_message("Writing Robinson-Foulds distances to " + rf_file)
        with open(rf_file, "w") as fout:
            fout.write("metric\tn_genes\trf_distance\n")
            for _, row in results:
                fout.write(row["metric"] + "\t" + str(row["n_genes"]) + "\t" + row["rf"] + "\n")

    print_message("Sensitivity analysis complete! Results in " + output_dir)


def run_metric(args):
    metric, gene_values, aln_dir, output_dir, n_top, model, ref_tree, iqtree_threads = args
    logs = []
    result = {"metric": metric, "n_genes": n_top, "rf": "NA"}

    # Sort: ascending for metrics where lower = better, descending otherwise
    ascending = metric in LOWER_IS_BETTER
    sorted_genes = sorted(gene_values, key=lambda x: x[1], reverse=not ascending)
    top_genes = [gene for gene, _ in sorted_genes[:n_top]]

    logs.append(print_message_str(
        metric + ": selected " + str(n_top) + " / " + str(len(gene_values)) + " genes"
    ))

    # Write gene list file (one alignment path per line)
    gene_list_file = os.path.join(output_dir, metric + ".gene_list.txt")
    with open(gene_list_file, "w") as fout:
        for gene in top_genes:
            fout.write(os.path.join(aln_dir, gene + ".trimmed.aln.fasta") + "\n")

    # Concatenate alignments with PhyKIT
    prefix = os.path.join(output_dir, metric)
    ret = subprocess.run(
        ["phykit", "create_concat", "-a", gene_list_file, "-p", prefix],
        capture_output=True, text=True
    )
    concat_fasta = prefix + ".fa"
    if ret.returncode != 0 or not os.path.isfile(concat_fasta):
        logs.append(print_message_str("Warning: phykit create_concat failed for metric " + metric))
        return logs, result

    logs.append(print_message_str(metric + ": alignment concatenated to " + concat_fasta))

    # Infer phylogenetic tree with IQ-TREE
    os.system(
        "iqtree3 --quiet -s " + concat_fasta +
        " -m " + model +
        " --fast -T " + str(iqtree_threads) + " >/dev/null 2>&1"
    )

    treefile = concat_fasta + ".treefile"
    if not os.path.isfile(treefile):
        logs.append(print_message_str("Warning: IQ-TREE did not produce a treefile for metric " + metric))
        return logs, result

    logs.append(print_message_str(metric + ": tree inferred at " + treefile))

    # Robinson-Foulds distance against reference tree
    if ref_tree:
        out = run_command(["phykit", "rf_dist", ref_tree, treefile])
        if out:
            result["rf"] = out.strip().split()[0]
            logs.append(print_message_str(
                metric + ": RF distance to reference = " + result["rf"]
            ))
        else:
            logs.append(print_message_str(
                "Warning: could not compute RF distance for metric " + metric
            ))

    result["n_genes"] = n_top
    return logs, result


def run_command(cmd):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return None
        return result.stdout
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def print_message_str(*message):
    return strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message))


def print_message(*message):
    print(print_message_str(*message))


if __name__ == "__main__":
    main()
