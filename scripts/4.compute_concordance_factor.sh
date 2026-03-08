#!/bin/bash
set -eou pipefail

#==============================================================================
# Title: 4.compute_concordance_factor.sh
# Description: This script computes concordance factors (gCF and sCF) 
#   using IQ-TREE based on 
#       - a supermatrix tree (SUPERMATRIX.treefile)
#       - a supermatrix alignment (SUPERMATRIX.aln.fasta)
#       - and a set of individual gene trees (ALL.trees).
#==============================================================================

IN_PHYLO_DIR="Phylogenomics"
SUPERMATRIX_ALN="SUPERMATRIX.aln.fasta"
SUPERMATRIX_TREE="SUPERMATRIX.treefile"
INDIVIDUAL_TREES="ALL.trees"
THREADS=8

conda activate busco

cd $IN_PHYLO_DIR

# Concordance factor compute
iqtree -t ${SUPERMATRIX_TREE} --gcf ${INDIVIDUAL_TREES} -s ${SUPERMATRIX_ALN} --scf 1000 -T ${THREADS}

conda deactivate