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
SUPERMATRIX_TREE="SUPERMATRIX.aln.fasta.treefile"
INDIVIDUAL_TREES="ALL.trees"
EVOL_MODEL="LG"
THREADS=64

conda activate busco-phylo

cd $IN_PHYLO_DIR

# Gene concordance factor (gCF)
# The --gcf option computes gene concordance factors based on the provided individual gene trees.
iqtree3 -t ${SUPERMATRIX_TREE} --gcf ${INDIVIDUAL_TREES}  -T 1 --prefix gcf

# Site concordance factor (sCF)
# The --scfl option computes site concordance factors based on the provided supermatrix alignment and tree,
#  using the specified evolutionary model.
iqtree3 -te ${SUPERMATRIX_TREE} -s ${SUPERMATRIX_ALN} --mset ${EVOL_MODEL} --scfl 1000 -T ${THREADS} --prefix scf

conda deactivate