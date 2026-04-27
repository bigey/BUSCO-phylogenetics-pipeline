#!/bin/bash
set -euo pipefail

#==============================================================================
# Set input and output directories, BUSCO parameters, and other settings.
# Modify these variables as needed for your specific use case.
#==============================================================================

# Directory containing genome assemblies in FASTA format (e.g.,SPECIES.fasta)
GENOME_DIR=Genomes

# Directory to store BUSCO results. 
# Will contain one subdirectory per genome
BUSCO_RESULTS="Busco-results"

# Choose a lineage dataset for BUSCO.
# See documentation for available lineages (https://busco.ezlab.org/busco_userguide.html)
#   or using the command line: busco --list-datasets
LINEAGE="saccharomycetes"

# Choose a gene predictor. 
# Options are: miniprot (default), metaeuk, augustus. 
#   miniprot is the fastest but may be less accurate, 
#   metaeuk is a good compromise between speed and accuracy,
#   augustus is the slowest but most accurate
# PREDICTOR="miniprot"
PREDICTOR="metaeuk"
# PREDICTOR="augustus"

# Number of threads to use for BUSCO. 
# Adjust based on your system's resources.
THREADS=16

#==============================================================================
# You don't need to modify anything below this line unless you want to change 
#   the BUSCO command or the way results are organized.
# 
# The script will run BUSCO for each genome sequences in `GENOME_DIR` 
#   and create links in the BUSCO result directories: `BUSCO_RESULTS`. 
#==============================================================================

mkdir -p ${BUSCO_RESULTS}

# Use a conda environment where all the required software are installed
conda activate busco-phylo

for genome in $GENOME_DIR/*.fasta; do
    PREFIX=$(basename ${genome} .fasta)
    OUT_DIR=busco_${PREFIX}

    echo "Processing genome ${genome}..."
    busco \
        --in ${genome} \
        --lineage_dataset ${LINEAGE} \
        --mode genome \
        --${PREDICTOR} \
        --out ${OUT_DIR} \
        --force \
        --cpu ${THREADS}
    echo "Finished processing genome ${genome}"

    # Clean up unnecessary directories and files
    find ${OUT_DIR} -name "logs" -o -name "tmp" -o -name "blast_db" | xargs rm -rf
    
    # Create a corresponding link to the BUSCO results directory
    cd ${BUSCO_RESULTS}
    ln -s ../${OUT_DIR}/run_${LINEAGE}_odb* run_${PREFIX}
    cd ..
done

conda deactivate
echo "done"

