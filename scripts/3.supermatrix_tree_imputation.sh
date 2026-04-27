#!/bin/bash
set -eou pipefail

#==============================================================================
# Title: 3.supermatrix_tree_imputation.sh
# Description: This script performs maximum likelihood (ML) phylogenetic tree 
#   imputation using IQ-TREE on a supermatrix alignment (SUPERMATRIX.aln.fasta) 
#   and generates a supermatrix tree (SUPERMATRIX.treefile)
#==============================================================================

PHYLO_DIR="Phylogenomics" # Directory containing the supermatrix alignment and where the tree will be generated
SUPERMATRIX_ALN="SUPERMATRIX.aln.fasta" # Supermatrix alignment file in FASTA format
PREFIX="SUPERMATRIX" # Prefix for output files (e.g., SUPERMATRIX.treefile, SUPERMATRIX.log, etc.)
EVOL_MODEL="LG"  # Evolutionary model for tree inference (e.g., LG, JTT, WAG, etc.)
THREADS=128 # Number of threads to use for IQ-TREE (adjust based on your system's capabilities)

conda activate busco-phylo

cd $PHYLO_DIR

# ML phylogenetic tree imputation
iqtree3 -s ${SUPERMATRIX_ALN} --prefix ${PREFIX} --mset ${EVOL_MODEL} --ufboot 1000 -T AUTO --threads-max ${THREADS}


conda deactivate