#!/usr/bin/env python3
import json
import logging
import sys
import argparse
import colorlog
from pathlib import Path


def extract_busco_fields(json_file):
    """Extract BUSCO fields from a single JSON file."""
    try:
        with open(json_file) as f:
            data = json.load(f)

        results = data.get("results", {})
        out = data.get("parameters", {}).get("out", "")

        return {
            "out": out,
            "Complete percentage": results.get("Complete percentage", ""),
            "Single pct": results.get("Single copy percentage", ""),
            "Multi pct": results.get("Multi copy percentage", ""),
            "Fragmented pct": results.get("Fragmented percentage", ""),
            "Missing pct": results.get("Missing percentage", ""),
        }
    except Exception as e:
        logger.warning(f"Failed to parse {json_file}: {e}")
        return None


def process_files(file_paths):
    """Process individual JSON files."""
    results = []
    for file_path in file_paths:
        path = Path(file_path)
        if not path.exists():
            logger.warning(f"File not found: {file_path}")
            continue
        if not path.suffix == ".json":
            logger.warning(f"Skipping non-JSON file: {file_path}")
            continue

        extracted = extract_busco_fields(path)
        if extracted:
            results.append(extracted)

    return results


def process_directory(directory, recursive=True):
    """Scan directory for BUSCO JSON files."""
    dir_path = Path(directory)
    if not dir_path.is_dir():
        logger.critical(f"Directory not found: {directory}")
        sys.exit(1)

    results = []
    pattern = "**/short_summary.specific.*.json" if recursive else "short_summary.specific.*.json"

    for json_file in dir_path.glob(pattern):
        extracted = extract_busco_fields(json_file)
        if extracted:
            results.append(extracted)

    return results


def print_results(results):
    """Print results as tab-separated values."""
    if not results:
        logger.warning("No BUSCO data extracted")
        return

    headers = ["Species", "Complete percentage", "Single pct", "Multi pct", "Fragmented pct", "Missing pct"]
    print("\t".join(headers))

    for result in results:
        values = [str(result.get("out" if h == "Species" else h, "")) for h in headers]
        print("\t".join(values))


def main(args):
    results = []

    if args.files:
        results.extend(process_files(args.files))

    if args.directory:
        for directory in args.directory:
            results.extend(process_directory(directory, recursive=not args.no_recursive))

    if not args.files and not args.directory:
        logger.critical("No input specified. Use --files or --directory.")
        sys.exit(1)

    print_results(results)


if __name__ == "__main__":
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
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler = colorlog.StreamHandler()
    handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.addHandler(handler)

    parser = argparse.ArgumentParser(
        description="Extract fields from BUSCO JSON files and output as tab-separated values."
    )
    parser.add_argument(
        "--files",
        nargs="+",
        help="BUSCO JSON files to process",
    )
    parser.add_argument(
        "--directory",
        nargs="+",
        help="Directories to scan for BUSCO JSON files",
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Don't scan subdirectories when processing directories",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )

    args = parser.parse_args()
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    main(args)
