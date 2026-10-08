#!/usr/bin/env python3

# Compute 7 PhyKIT information-content metrics for each single-copy BUSCO gene.
# Metrics: alignment length, average bipartition support, relative composition
# variability, median long-branch score, saturation, treeness, treeness/RCV.
# Output is a tab-delimited file with columns: gene, metric, value.
#
# Step 1 of the BUSCO phylogenomics sensitivity analysis pipeline.
# Reference: https://jlsteenwyk.com/tutorials/ub_sensitivity_analysis.html#Section2
#
# Dependencies:
#   - PhyKIT

import os
import sys
import time
import shutil
import logging
import colorlog
import argparse
import subprocess
import multiprocessing as mp


# External dependencies
external_dependencies = ["phykit"]


def check_dependency(name):
    if shutil.which(name) is None:
        logger.critical(f"Program {name} is required but not in the PATH!")
        logger.critical(f"Please install {name}.")
        sys.exit(1)


def main(args):
    aln_dir = os.path.abspath(args.trimmed_alignments)
    tree_dir = os.path.abspath(args.trees)
    output_file = os.path.abspath(args.output)
    threads = args.threads

    # Validate input directories
    if not os.path.isdir(aln_dir):
        logger.critical(f"{aln_dir} is not a directory!")
        sys.exit(1)
    if not os.path.isdir(tree_dir):
        logger.critical(f"{tree_dir} is not a directory!")
        sys.exit(1)

    # Discover genes: pair each trimmed alignment with its gene tree
    logger.info(f"Scanning {aln_dir} for trimmed alignments...")
    jobs = []
    missing_trees = []

    for fname in sorted(os.listdir(aln_dir)):
        if not fname.endswith(".trimmed.aln.fasta"):
            continue
        busco_id = fname.replace(".trimmed.aln.fasta", "")
        aln_path = os.path.join(aln_dir, fname)
        tree_path = os.path.join(tree_dir, fname + ".treefile")

        if not os.path.isfile(tree_path):
            missing_trees.append(busco_id)
            continue

        jobs.append((busco_id, aln_path, tree_path))

    if missing_trees:
        logger.warning(f"{len(missing_trees)} genes have no treefile and will be skipped")

    if not jobs:
        logger.critical("No gene/tree pairs found. Exiting.")
        sys.exit(1)

    logger.info(f"{len(jobs)} genes found with paired trimmed alignment and gene tree")

    # Compute metrics in parallel
    logger.info(f"Computing 7 PhyKIT metrics per gene using {threads} threads...")
    pool = mp.Pool(processes=threads)
    results = pool.map(compute_metrics, jobs)
    pool.close()
    pool.join()

    # Write results
    logger.info(f"Writing results to {output_file}")
    n_ok = 0
    with open(output_file, "w") as fout:
        for entry in results:
            if entry is None:
                continue
            busco_id, metrics = entry
            for metric, value in metrics.items():
                fout.write(busco_id + "\t" + metric + "\t" + value + "\n")
            n_ok += 1

    logger.info(f"{n_ok} genes written to {output_file}")
    logger.info("Step 1 complete!")


def compute_metrics(args):
    busco_id, aln_path, tree_path = args
    metrics = {}

    # 1. Alignment length
    out = run_command(["phykit", "aln_len", aln_path])
    metrics["aln_len"] = out.strip().split()[0] if out else "NA"

    # 2. Average bipartition support (mean from multi-line output)
    out = run_command(["phykit", "bss", tree_path])
    metrics["abs"] = "NA"
    if out:
        for line in out.splitlines():
            if "mean" in line.lower():
                metrics["abs"] = line.split()[-1]
                break

    # 3. Relative composition variability
    out = run_command(["phykit", "rcv", aln_path])
    metrics["rcv"] = out.strip().split()[0] if out else "NA"

    # 4. Median long-branch score (median from multi-line output)
    out = run_command(["phykit", "lbs", tree_path])
    metrics["lbs"] = "NA"
    if out:
        for line in out.splitlines():
            if "median" in line.lower():
                metrics["lbs"] = line.split()[-1]
                break

    # 5. Treeness
    out = run_command(["phykit", "tness", tree_path])
    metrics["treeness"] = out.strip().split()[0] if out else "NA"

    # 6. Saturation (second field: absolute value of saturation minus 1)
    out = run_command(["phykit", "sat", "-a", aln_path, "-t", tree_path])
    metrics["saturation"] = "NA"
    if out:
        parts = out.strip().split()
        if len(parts) >= 2:
            metrics["saturation"] = parts[1]

    # 7. Treeness/RCV ratio (first field)
    out = run_command(["phykit", "tor", "-a", aln_path, "-t", tree_path])
    metrics["treeness_over_rcv"] = "NA"
    if out:
        parts = out.strip().split()
        if parts:
            metrics["treeness_over_rcv"] = parts[0]

    return busco_id, metrics


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
        description="Compute PhyKIT information-content metrics for each single-copy BUSCO gene"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Turn on verbose mode."
    )
    parser.add_argument(
        "--trimmed-alignments", type=str, required=True,
        help="Directory containing trimmed alignment files (*.trimmed.aln.fasta)"
    )
    parser.add_argument(
        "--trees", type=str, required=True,
        help="Directory containing individual gene trees (*.trimmed.aln.fasta.treefile)"
    )
    parser.add_argument(
        "--output", type=str, default="info_content_genes.txt",
        help="Output tab-delimited file (default: info_content_genes.txt)"
    )
    parser.add_argument(
        "--threads", type=int, default=8,
        help="Number of parallel threads (default: 8)"
    )
    args = parser.parse_args()

    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    # handler configuration
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

    # logger configuration
    logger = logging.getLogger()
    logger.addHandler(handler)
    logger.setLevel(log_level)

    start_time = time.time()

    # Check dependencies
    list(map(check_dependency, external_dependencies))

    main(args)

    end_time = time.time()
    execution_time = time.strftime("%Hh:%Mm:%Ss", time.gmtime(end_time - start_time))
    logger.info(f"Execution time: {execution_time}")
    logger.info("Done")
