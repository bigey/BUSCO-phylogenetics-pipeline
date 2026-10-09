#!/usr/bin/env python3

# Utility script to construct species phylogenies using single-copy BUSCO proteins.
# Can perform ML supermatrix phylogeny or generate datasets for supertree methods.
# Works directly from BUSCO output, as long as the same BUSCO dataset
# has been used for each genome

# Python modules
import argparse
import logging
import colorlog
import multiprocessing as mp
import subprocess
import shutil
import time
import os
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# External dependency
external_dependencies = ["mafft","clipkit","iqtree3"]

def check_dependency(name):
    if shutil.which(name) is None:
        logger.critical(f"Program {name} is required bu not in the PATH!")
        logger.critical(f"Please install {name}.")
        sys.exit(1)

def main(args):

    start_directory = os.path.abspath(args.busco_dir)
    working_directory = os.path.abspath(args.output_dir)
    threads = int(args.threads)
    supermatrix = args.supermatrix
    supertree = args.supertree
    concordance = args.concordance
    outgroup = args.outgroup
    model = args.model
    percent_single_copy = args.percent_single_copy
    stop_early = args.stop_early

    logger.info("Checking input parameters...")

    # Check input parameters
    if stop_early and concordance:
        print("Error! '--concordance' cannot be used with '--stop-early'")
        sys.exit(1)
    elif stop_early and not supermatrix:
        print("Error! '--stop-early' can only be used with '--supermatrix'")
        sys.exit(1)
    elif stop_early and supertree:
        print("Error! '--stop-early' cannot be used with '--supertree'")
        sys.exit(1)

    if concordance:
        supermatrix = True
        supertree = True
        logger.info("Option '--concordance' required automatic selection of both '--supermatrix' and '--supertree' options.")

    if not supermatrix and not supertree:
        print("Error! Please select at least one of '--supermatrix' or '--supertree'")
        sys.exit(1)

    # Check input directory exists
    if not os.path.isdir(start_directory):
        print("Error! " + start_directory + " is not a directory!")
        sys.exit(1)

    # Check if output directory already exists
    if os.path.isdir(working_directory):
        logger.warning(f"Directory {working_directory} already exists!")
        os.chdir(working_directory)

        if not (os.path.isdir("proteins") and os.path.isdir("alignments") and os.path.isdir("trimmed_alignments")):
            logger.critical("One required output directory is missing!")
            logger.critical("See missing directory in: proteins/, alignments/, trimmed_alignments/, trees/")
            sys.exit(1)
    else:
        os.makedirs(working_directory)
        os.chdir(working_directory)
        os.makedirs("proteins", exist_ok=True)
        os.makedirs("alignments", exist_ok=True)
        os.makedirs("trimmed_alignments", exist_ok=True)
        os.makedirs("trees", exist_ok=True)

    # Starting pipeline
    logger.info("Starting phylogenomics pipeline...")

    # Align and trim BUSO single copy proteins
    os.chdir(working_directory)

    if not os.path.isfile("busco.aligned-files.txt"):
        logger.info("Starting constructing a supermatrix file...")

        # Scan directory to identify BUSCO runs (directories should begin with 'run_')
        os.chdir(start_directory)
        busco_dirs = []

        for item in os.listdir("."):
            if item[0:4] == "run_":
                if os.path.isdir(item):
                    busco_dirs.append(item)

        nb_run = str(len(busco_dirs))
        logger.info(f"{nb_run} BUSCO runs were found in {start_directory}")

        buscos = {}
        all_species = []

        # Parse BUSCO sequences from each run
        for directory in busco_dirs:
            os.chdir(start_directory)

            species = directory.split("run_")[1]
            all_species.append(species)

            os.chdir(directory)
            os.chdir("busco_sequences")
            os.chdir("single_copy_busco_sequences")

            for busco in os.listdir("."):
                if busco.endswith(".faa"):
                    busco_name = busco[0 : len(busco) - 4]
                    record = SeqIO.read(busco, "fasta")
                    new_record = SeqRecord(Seq(str(record.seq)), id=species, description="")

                    if busco_name not in buscos:
                        buscos[busco_name] = []

                    buscos[busco_name].append(new_record)

        nb_single_copy = str(len(buscos))
        logger.info(f"{nb_single_copy} BUSCO single-copy genes were found")

        # Test if the species given as outgroup are in the "all_species" list
        if outgroup is not None:
            outgroup = outgroup.replace(" ", "")
            for sp in outgroup.split(","):
                if sp not in all_species:
                    logger.critical(f"Outgroup species {sp} is not in the input species list")
                    logger.critical(f"List: {all_species}")
                    sys.exit(1)

        single_copy_buscos = []

        # Select BUSCO genes that are present in all species
        if percent_single_copy == 1.0:
            logger.info(f"Identifying BUSCO genes that are present in all species...")

            for busco in buscos:
                if len(buscos[busco]) == len(all_species):
                    single_copy_buscos.append(busco)

            if len(single_copy_buscos) == 0:
                logger.critical("No BUSCO genes were present in all species! Exiting")
                sys.exit(0)
            else:
                nb_single_copy_gene = str(len(single_copy_buscos))
                nb_species = str(len(all_species))
                logger.info(f"{nb_single_copy_gene} BUSCO genes are present in all {nb_species}.")
        
        # Select BUSCO genes that are present in at least a given faction of species
        else:
            for busco in buscos:
                percent_species_with_single_copy = len(buscos[busco]) / (len(all_species) * 1.0)

                if percent_species_with_single_copy >= percent_single_copy:
                    single_copy_buscos.append(busco)
            
            nb_single_copy_gene = len(single_copy_buscos)
            logger.info(f"{nb_single_copy_gene} BUSCO genes are present in at least {percent_single_copy} of species")

        # Write BUSCO sequences to protein folder
        logger.info(f"Writing BUSCO protein sequences to: {os.path.join(working_directory, 'proteins')}")

        for busco in single_copy_buscos:
            busco_seqs = buscos[busco]
            SeqIO.write(
                busco_seqs,
                os.path.join(working_directory, "proteins", busco + ".faa"),
                "fasta",
            )

        # Align each BUSCO family
        logger.info(f"Aligning protein sequences using {threads} threads, results to: {os.path.join(working_directory, 'alignments')}")
        mp_commands = []

        for busco in single_copy_buscos:
            mp_commands.append(
                [
                    os.path.join(working_directory, "proteins", busco + ".faa"),
                    os.path.join(working_directory, "alignments", busco + ".aln.fasta"),
                ]
            )

        pool = mp.Pool(processes=threads)
        pool.map(run_mafft, mp_commands)

        logger.info("All alignment jobs finished!")

        # Trim resulting alignments
        logger.info(f"Trimming alignments using {threads} threads, results to: {os.path.join(working_directory, 'trimmed_alignments')}")
        mp_commands = []

        for busco in single_copy_buscos:
            mp_commands.append(
                [
                    os.path.join(working_directory, "alignments", busco + ".aln.fasta"),
                    os.path.join(
                        working_directory,
                        "trimmed_alignments",
                        busco + ".trimmed.aln.fasta",
                    ),
                ]
            )

        pool = mp.Pool(processes=threads)
        pool.map(run_clipkit, mp_commands)

        # Save the aligned sequence file path
        os.chdir(working_directory)
        trim_files = []

        for busco in single_copy_buscos:
            trim_files.append(os.path.join(working_directory, "trimmed_alignments", busco + ".trimmed.aln.fasta"))

        with open("busco.aligned-files.txt", "w") as f_path:   
            f_path.write("\n".join(trim_files) + "\n")

        logger.info("All trimming jobs finished!")

    else:
        logger.warning("All BUSCO proteins are already aligned and trimmed! Skipping this step...")
        os.chdir(working_directory)

        all_species = []
        single_copy_buscos = []
        all_species_set = set()

        with open("busco.aligned-files.txt") as f_path:

            for line in f_path:
                trimmed_path = line.strip()

                if not trimmed_path:
                    continue
                if not os.path.isfile(trimmed_path):
                    logger.critical(f"Expected trimmed alignment file not found: {trimmed_path}")
                    sys.exit(1)
                
                filename = os.path.basename(trimmed_path)

                if filename.endswith(".trimmed.aln.fasta"):
                    busco_name = filename[: -len(".trimmed.aln.fasta")]
                elif filename.endswith(".aln.fasta"):
                    busco_name = filename[: -len(".aln.fasta")]
                else:
                    busco_name = os.path.splitext(filename)[0]

                single_copy_buscos.append(busco_name)

                for record in SeqIO.parse(trimmed_path, "fasta"):
                    all_species_set.add(str(record.id))

        all_species = sorted(all_species_set)
    
    # Compute supermatrix tree phylogeny
    if supermatrix:
        logger.info("Creation of a supermatrix species tree was selected.")
        os.chdir(working_directory)

        # Test if supermatrix tree file already exist
        if not os.path.isfile("SUPERMATRIX.aln.fasta.treefile"):
            logger.info("Starting supermatrix analysis...")
            logger.info("Concatenating all trimmed alignments...")

            os.chdir(os.path.join(working_directory, "trimmed_alignments"))
            alignments = {}

            # Initialize alignments dictionary with all species as keys and empty strings as values
            for species in all_species:
                alignments[species] = ""

            # If percent_single_copy is 1.0, we can simple just concatenate alignments
            if percent_single_copy == 1.0:
                for alignment in os.listdir("."):
                    # For each alignment file,
                    #   append the sequence to the corresponding species 
                    #   in the alignments dictionary
                    for record in SeqIO.parse(alignment, "fasta"):
                        alignments[str(record.id)] += str(record.seq)

            # Else, we need to check if a species is missing from a family,
            # if so append with "-" to represent missing data
            else:
                for alignment in os.listdir("."):
                    # Keep track of which species are present or missing
                    missing_species = all_species[:]

                    for record in SeqIO.parse(alignment, "fasta"):
                        alignments[str(record.id)] += str(record.seq)
                        missing_species.remove(str(record.id))

                    if len(missing_species) > 0:
                        # There are missing species,
                        #   we need to fill in the alignment with "-"
                        #   for those species to represent missing data

                        # Get the length of the first alignment
                        first_seq = SeqIO.parse(alignment, "fasta")
                        seq_len = len(str(first_seq.seq))

                        for species in missing_species:
                            #  fill with "-" character
                            alignments[species] += "-" * seq_len

            # Write supermatrix alignment to file
            logger.info(f"Writing supermatrix file to: {os.path.join(working_directory, 'SUPERMATRIX.aln.fasta')}")
            os.chdir(working_directory)
            seq_records = []

            for species in alignments:
                seq_records.append(SeqRecord(Seq(alignments[species]), id = species))
            
            SeqIO.write(seq_records, "SUPERMATRIX.aln.fasta", "fasta")

            # All alignments should be the same length, 
            #   so we can just check the length of the first one
            logger.info(f"Supermatrix alignment is {len(seq_records[0].seq)} amino acids in length")

            if stop_early:
                logger.info("Stopping here as requested with --stop-early option. Exiting.")
                sys.exit(0)

            # Imputing phylogenomic tree
            logger.info(f"Start species tree imputation using {threads} threads...")
            logger.info(f"Species tree will go to: {os.path.join(working_directory, 'SUPERMATRIX.aln.fasta.treefile')}")

            # Test if outgroup option is set
            if outgroup is None:
                subprocess.call(["iqtree3", "--undo", "-s", "SUPERMATRIX.aln.fasta", "--quiet", "-m", model, "-T", str(threads), "--ufboot", "1000"])
            else:
                subprocess.call(["iqtree3", "--undo", "-s", "SUPERMATRIX.aln.fasta", "--quiet", "-o", outgroup, "-m", model, "-T", str(threads), "--ufboot", "1000"])

            logger.info("Supermatrix species tree imputation finished!")

        else:
            logger.warning("Supermatrix tree already present. Skiping this step...")

    # Compute supertree phylogeny
    if supertree:
        logger.info("Creation of a supertree was selected.")
        os.chdir(working_directory)

        # Test if concatenation tree file already exist
        if not os.path.isfile(os.path.join(working_directory, "ALL.trees")):

            logger.info("Starting supertree analysis...")
            logger.info(f"Generating a phylogenetic tree for each BUSCO gene/protein, {threads} threads used...")
            iqtree_commands = []

            for busco in single_copy_buscos:
                iqtree_commands.append(
                    [
                        os.path.join("trimmed_alignments", busco + ".trimmed.aln.fasta"),
                        model, outgroup,
                    ]
                )

            pool = mp.Pool(processes=threads)
            pool.map(run_iqtree, iqtree_commands)

            # Move all IQ-TREE generated files to trees folder
            os.makedirs("trees/iqtree_files", exist_ok=True)

            for f in Path("trimmed_alignments").glob("*.treefile"):
                shutil.move(str(f), "trees/")

            for f in Path("trimmed_alignments").glob("*.trimmed.aln.fasta.*"):
                shutil.move(str(f), "trees/iqtree_files")

            logger.info("All tree imputation jobs are finished")
            logger.info(f"Concatenating all trees to: {os.path.join(working_directory, 'ALL.trees')}")

            with open("ALL.trees", "w") as out_f:
                for tree_path in sorted(Path("trees").glob("*.treefile")):
                    with open(tree_path) as in_f:
                        out_f.write(in_f.read())
        else:
            logger.warning("Concatenation tree file already exists. Skiping this step...")
    
    if concordance:
        logger.info("Computing gene/site concordance factors was selected.")

        if not os.path.isfile(os.path.join(working_directory, "concordance-factors.cf.tree")):
            os.chdir(working_directory)

            # Calculate gCF and sCF using IQ-TREE
            subprocess.call(["iqtree3", "--undo", "--quiet", "-te", "SUPERMATRIX.aln.fasta.treefile", "-s", "SUPERMATRIX.aln.fasta", "--gcf", "ALL.trees", "-m", model, "--scf", "1000", "--prefix", "concordance-factors"])

            logger.info("gCF and sCF estimation complete.")
            logger.info(f"See annotated tree file: {os.path.join(working_directory, 'concordance-factors.cf.tree')}")
        else:
            logger.warning("Concordance factor annotated tree already exists. Skiping this step...")

    # Final message
    logger.info("BUSCO phylogenomics pipeline complete!")

def run_mafft(io):
    with open(io[1], "w") as fout:
        subprocess.call(["mafft", "--quiet", "--thread", "1", io[0]], stdout=fout)

def run_clipkit(io):
    subprocess.call(["clipkit", io[0], "--quiet", "--mode", "smart-gap", "--output", io[1]])

def run_iqtree(param):
    alignment = param[0]
    model = param[1]
    outgrp = param[2]

    if outgrp is None:
        subprocess.call(["iqtree3", "--undo", "--quiet", "-s", alignment, "-m", model, "-T", "1"])
    else:
        subprocess.call(["iqtree3", "--undo", "--quiet", "-s", alignment, "-m", model, "-o", outgrp, "-T", "1"])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Perform phylogenomic reconstruction using single-copy BUSCO genes."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Turn on verbose mode."
    )
    parser.add_argument(
        "--supermatrix",
        action="store_true",
        help="Concatenate alignments of single-copy BUSCO proteins and perform supermatrix ML phylogeny using IQ-TREE",
    )
    parser.add_argument(
        "--supertree",
        action="store_true",
        help="Generate individual ML phylogenies of each BUSCO proteins using IQ-TREE. Then concatenate all trees into a supertree file",
    )
    parser.add_argument(
        "--concordance",
        action="store_true",
        default=False,
        help="Calculate concordance factors (gCF and sCF) for phylogenetic trees, automatic turn on --supermatrix and --supertree if not already selected",
    )
    parser.add_argument(
        "--stop-early",
        action="store_true",
        default=False,
        help="Stop pipeline early after generating datasets (before phylogeny inference), only relevant if --supermatrix method is selected, incompatible with --concordance",
    )
    parser.add_argument(
        "--percent-single-copy",
        type=float,
        default=1.0,
        help="Only keep BUSCO genes present in single copy in at least this fraction of the species (default 1.0, i.e. 100%)",
    )
    parser.add_argument(
        "--outgroup",
        type=str,
        default=None,
        help="List of outgroup species separated by a comma (default None)",
        required=False,
    )
    parser.add_argument(
        "--model",
        type=str,
        default="LG+R4+F",
        help="Protein evolution model to use for tree inference (see IQ-TREE documentation, default to LG+R4+F model)",
    )
    parser.add_argument(
        "--threads", type=int, default=8, help="Number of threads to use (default 8)"
    )
    parser.add_argument(
        "--busco-dir",
        type=str,
        help="Directory containing completed BUSCO runs",
        required=True,
    )
    parser.add_argument(
        "--output-dir", type=str, help="Output directory to store results", required=True
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
