#!/bin/bash
set -euo pipefail

REFERENCE_TREE="Phylogenomics/SUPERMATRIX.aln.fasta.treefile"
TRIMMED_ALIGNMENTS="Phylogenomics/trimmed_alignments"
TREES_DIR="Phylogenomics/trees"
OUTPUT_DIR="Phylogenomics"
MODEL="LG+R4+F"
TOP_FRACTION="0.5"
THREADS=64

conda activate busco-phylo

# Run the script to compute the information content of each gene tree
# This script will read the trimmed alignments and their corresponding gene trees, 
#    compute several metrics related to the information content of each gene, 
#    such as the number of informative sites, tree length, and other relevant statistics. 
# The results will be saved in a TSV file for further analysis.
python3 compute-gene-metrics.py \
    --trimmed-alignments ${TRIMMED_ALIGNMENTS} \
    --trees ${TREES_DIR} \
    --output ${OUTPUT_DIR}/gene-metrics.tsv \
    --threads ${THREADS}

# We generate a PDF file with metric distributions.
python3 plot-gene-metrics.py --input ${OUTPUT_DIR}/gene-metrics.tsv --output ${OUTPUT_DIR}/gene-metrics.pdf

# Determine if the phylogenetic trees inferred using subsampling differ from the species/referece tree. 
# Will perform a sensitivity analysis to assess how the choice of genes 
#    (based on their information content) affects the inferred phylogenetic tree. 
# Select the top 0.X fraction of best scoring genes based on their information content metrics,
#    and infer a new phylogenetic tree using only these genes.
# Then compare this tree to the reference tree by using the Robinson-Foulds metric, 
#    which quantifies the distance between two trees.
#    Robinson & Foulds, Mathematical Biosciences (1981), doi: 10.1016/0025-5564(81)90043-2.
python3 gene-sensitivity-analysis.py \
    --ref-tree ${REFERENCE_TREE} \
    --input ${OUTPUT_DIR}/gene-metrics.tsv \
    --trimmed-alignments ${TRIMMED_ALIGNMENTS} \
    --output-dir ${OUTPUT_DIR}/gene-sensitivity-top-${TOP_FRACTION} \
    --top-fraction ${TOP_FRACTION} \
    --model ${MODEL} \
    --threads ${THREADS}

conda deactivate
