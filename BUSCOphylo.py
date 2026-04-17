#!/usr/bin/env python3

# Utility script to construct species phylogenies using single-copy BUSCO proteins.
# Can perform ML supermatrix phylogeny or generate datasets for supertree methods.
# Works directly from BUSCO output, as long as the same BUSCO dataset
# has been used for each genome
#
# Dependencies:
#   - BioPython
#   - MAFFT
#   - ClipKit
#   - IQ-TREE

# from operator import mod
import os
import sys
import argparse
import multiprocessing as mp

from time import gmtime, strftime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description="Perform phylogenomic reconstruction using single-copy BUSCO proteins")

    parser.add_argument("--supermatrix", action="store_true", help="Concatenate alignments of single-copy BUSCO proteins and perform supermatrix ML phylogeny using IQ-TREE")
    parser.add_argument("--supertree", action="store_true", help="Generate individual ML phylogenies of each BUSCO proteins using IQ-TREE. Then concatenate all trees into a supertree file")
    parser.add_argument("--concordance", action="store_true", default=False, help="Calculate concordance factors (gCF and sCF) for phylogenetic trees, automatic turn on --supermatrix and --supertree if not already selected")
    parser.add_argument("--stop-early", action="store_true", default=False, help="Stop pipeline early after generating datasets (before phylogeny inference), only relevant if --supermatrix method is selected, incompatible with --concordance")
    parser.add_argument("--percent-single-copy", type=float, default=1.0, help="Only keep BUSCO genes present in single copy in at least this fraction of the species (default 1.0, i.e. 100%)")
    parser.add_argument("--model", type=str, default="LG", help="Protein evolution model to use for IQ-TREE phylogeny inference (see IQ-TREE documentation, default to LG model)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use (default 8)")
    parser.add_argument("--directory", type=str, help="Directory containing completed BUSCO runs", required=True)
    parser.add_argument("--output", type=str, help="Output directory to store results", required=True)
    args = parser.parse_args()

    start_directory = os.path.abspath(args.directory)
    working_directory = os.path.abspath(args.output)
    threads = int(args.threads)
    supermatrix = args.supermatrix
    supertree = args.supertree
    concordance = args.concordance
    model = args.model
    percent_single_copy = args.percent_single_copy
    stop_early = args.stop_early
  
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
        print("Selection of the '--concordance' option required automatic selection of both '--supermatrix' and '--supertree' options.")

    if not supermatrix and not supertree:
        print("Error! Please select at least one of '--supermatrix' or '--supertree'")
        sys.exit(1)

    if model not in ["Model","Blosum62","cpREV","Dayhoff","DCMut","EAL","ELM","FLAVI","FLU","GTR20","HIVb","HIVw","JTT","JTTDCMut","LG","mtART","mtMAM","mtREV","mtZOA","mtMet","mtVer","mtInv","NQ.bird","NQ.insect","NQ.mammal","NQ.pfam","NQ.plant","NQ.yeast","Poisson","PMB","Q.bird","Q.insect","Q.mammal","Q.pfam","Q.plant","Q.yeast","rtREV","VT","WAG"]:
        print("Error! Specified model " + model + " is not supported by IQ-TREE. See IQ-TREE documentation for supported models.")
        sys.exit(1)
    
    # Check input directory exists
    if not os.path.isdir(start_directory):
        print("Error! " + start_directory + " is not a directory!")
        sys.exit(1)

    # Check if output directory already exists
    if os.path.isdir(working_directory):
        print("Error! " + working_directory + " already exists")
        sys.exit(1)
    else:
        os.mkdir(working_directory)

    # Begin processing BUSCO runs
    print_message("Starting phylogenomics pipeline...")

    # Scan directory to identify BUSCO runs (directories should begin with 'run_')
    os.chdir(start_directory)
    busco_dirs = []

    for item in os.listdir("."):
        if item[0:4] == "run_":
            if os.path.isdir(item):
                busco_dirs.append(item)

    print_message(str(len(busco_dirs)) + " BUSCO runs were found in " + start_directory)

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
                busco_name = busco[0:len(busco) - 4]
                record = SeqIO.read(busco, "fasta")
                new_record = SeqRecord(Seq(str(record.seq)), id=species, description="")

                if busco_name not in buscos:
                    buscos[busco_name] = []

                buscos[busco_name].append(new_record)

    print_message((str(len(buscos))) + " BUSCOs were found")

    # Select BUSCO genes that are present in all species
    single_copy_buscos = []

    if percent_single_copy == 1.0:
        print_message("Identifying BUSCOs that are single copy in all " + str(len(all_species)) + " species")

        for busco in buscos:
            if len(buscos[busco]) == len(all_species):
                single_copy_buscos.append(busco)

        if len(single_copy_buscos) == 0:
            print_message("0 BUSCO families were present and single copy in all species! Exiting")
            sys.exit(0)
        else:
            print_message(str(len(single_copy_buscos)) + " BUSCOs are single copy in all " + str(len(all_species)) + " species")
    else:
        # Identify BUSCOs that are single copy and present in above % of species
        psc = percent_single_copy

        for busco in buscos:
            percent_species_with_single_copy = (len(buscos[busco]) / (len(all_species) * 1.0))

            if percent_species_with_single_copy >= psc:
                single_copy_buscos.append(busco)

        print_message(str(len(single_copy_buscos)) + " BUSCOs are single copy in >= " + str(psc) + " of species")

    os.chdir(working_directory)
    os.mkdir("proteins")
    os.mkdir("alignments")
    os.mkdir("trimmed_alignments")
    os.mkdir("trees")

    # Write BUSCO sequences to protein folder
    print_message("Writing BUSCO protein sequences to: " + os.path.join(working_directory, "proteins"))
    for busco in single_copy_buscos:
        busco_seqs = buscos[busco]
        SeqIO.write(busco_seqs, os.path.join(working_directory, "proteins", busco + ".faa"), "fasta")

    # Align each BUSCO family
    print_message("Aligning protein sequences using MAFFT to: ", os.path.join(working_directory))
    mp_commands = []
    for busco in single_copy_buscos:
        mp_commands.append([os.path.join(working_directory, "proteins", busco + ".faa"), os.path.join(working_directory, "alignments", busco + ".aln.fasta")])

    pool = mp.Pool(processes=threads)
    pool.map(run_mafft, mp_commands)

    # Trim alignments with ClipKit
    print_message("Trimming alignments using ClipKit to: ", os.path.join(working_directory, "trimmed_alignments"))
    mp_commands = []
    for busco in single_copy_buscos:
        mp_commands.append([os.path.join(working_directory, "alignments", busco + ".aln.fasta"), os.path.join(working_directory, "trimmed_alignments", busco + ".trimmed.aln.fasta")])

    pool = mp.Pool(processes=threads)
    pool.map(run_clipkit, mp_commands)

    print_message("Alignment and trimming complete!")

    # Compute Supermatrix
    if supermatrix:
        print_message("Beginning SUPERMATRIX analysis...")
        print_message("Concatenating all trimmed alignments for SUPERMATRIX analysis...")

        os.chdir(os.path.join(working_directory, "trimmed_alignments"))
        alignments = {}

        # Initialize alignments dictionary with all species as keys and empty strings as values
        for species in all_species:
            alignments[species] = ""

        # If percent_single_copy is 1.0, we can simple just concatenate alignments
        if percent_single_copy ==1.0:
            for alignment in os.listdir("."):
                # For each alignment file, 
                #  append the sequence to the corresponding species in the alignments dictionary
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
                    #  we need to fill in the alignment with "-" 
                    #  for those species to represent missing data
                    
                    # Get the length of the alignment
                    first = next(SeqIO.parse(alignment, "fasta")).seq
                    seq_len = len(str(first))
                     
                    for species in missing_species:
                        #  Fill with N * "-" 
                        alignments[species] += ("-" * seq_len)

        # Write supermatrix alignment to file
        os.chdir(working_directory)
        fo = open("SUPERMATRIX.aln.fasta", "w")
        for species in alignments:
            fo.write(">" + species + "\n")
            fo.write(alignments[species] + "\n")
        fo.close()
        
        # All alignments should be the same length, so we can just check the length of the first one
        first_align = next(iter(alignments.values()))
        print_message("Supermatrix alignment is " + str(len(first_align)) + " amino acids in length")

        if stop_early:
            print_message("Stopping here as requested with --stop-early. Exiting.")
            sys.exit(0)

        print_message("Reconstructing species phylogeny using IQ-TREE: tree will go to SUPERMATRIX.aln.fasta.treefile")
        os.system(f"iqtree3 -s SUPERMATRIX.aln.fasta --quiet --mset {model} --ufboot 1000 -T AUTO --threads-max {threads} >/dev/null 2>&1")

        print_message("SUPERMATRIX phylogeny construction complete! See treefile: SUPERMATRIX.aln.fasta.treefile")

    # Compute Supertree
    if supertree:
        print_message("Beginning SUPERTREE analysis...")
        print_message("Generating phylogenies using IQ-TREE for each BUSCO family for SUPERTREE analysis...")

        mp_commands = []

        for busco in single_copy_buscos:
            mp_commands.append([os.path.join("trimmed_alignments", busco + ".trimmed.aln.fasta"), model])

        pool = mp.Pool(processes=threads)
        pool.map(run_iqtree, mp_commands)

        # Move all IQ-TREE generated files to trees folder
        os.chdir(working_directory)
        os.mkdir("trees/iqtree_files")
        os.system("mv trimmed_alignments/*.treefile trees")
        os.system("mv trimmed_alignments/*.trimmed.aln.fasta.* trees/iqtree_files")

        print_message("Concatenating all trees to: ", os.path.join(working_directory, "ALL.trees"))
        os.system("cat trees/*.treefile > ALL.trees")

        print_message("Finished generating all individual trees.")

    if concordance:
        print_message("Calculating gene and site concordance factors (gCF and sCF) for SUPERMATRIX tree...")
        os.chdir(working_directory)
        
        print_message("Calculating gCF...")
        os.system("iqtree3 --quiet -t  SUPERMATRIX.aln.fasta.treefile --gcf ALL.trees --prefix gcf -T 1 >/dev/null 2>&1")
        
        print_message("Calculating sCF...")
        os.system(f"iqtree3 --quiet -te SUPERMATRIX.aln.fasta.treefile -s SUPERMATRIX.aln.fasta --scfl 1000 --prefix scf -T AUTO --threads-max {threads} >/dev/null 2>&1")

        print_message("gCF and sCF calculation complete! See gcf.*/scf.* for annotated treefile.")

    # Final message
    print_message("BUSCO phylogenomics pipeline complete!")

def run_mafft(io):
    os.system("mafft --quiet --thread 1 " + io[0] + " > " + io[1])

def run_clipkit(io):
    os.system("clipkit " + io[0] + " --quiet --mode smart-gap --output " + io[1])

def run_iqtree(param):
    alignment = param[0]
    model = param[1]
    os.system(f"iqtree3 --quiet -s {alignment} --mset {model} --ufboot 1000 -T 1")

def print_message(*message):
    print(strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message)))

if __name__ == "__main__":
    main()
