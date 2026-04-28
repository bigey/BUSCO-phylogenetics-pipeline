#!/bin/bash
set -euo pipefail

TRIMMED_ALIGNMENTS="Phylogenomics/trimmed_alignments"
TREES_DIR="Phylogenomics/trees"
OUTPUT_DIR="Phylogenomics"
MODEL="LG+R4+F"
TOP_FRACTION="0.50"
THREADS=64

conda activate busco-phylo

# Run the script to compute the information content of each gene tree
# This script will read the trimmed alignments and their corresponding gene trees, 
#    compute several metrics related to the information content of each gene, 
#    such as the number of informative sites, tree length, and other relevant statistics. 
# The results will be saved in a TSV file for further analysis.
# python3 compute-gene-metrics.py \
#     --trimmed-alignments ${TRIMMED_ALIGNMENTS} \
#     --trees ${TREES_DIR} \
#     --output ${OUTPUT_DIR}/gene-metrics.tsv \
#     --threads ${THREADS}

# We plot the results and generate a PDF file with the visualizations of the metric distributions.
# python3 plot-gene-metrics.py --input ${OUTPUT_DIR}/gene-metrics.tsv --output ${OUTPUT_DIR}/gene-metrics.pdf

# Determine if the phylogenetic trees inferred using subsampling differ from the species/referece tree. 
# We will perform a sensitivity analysis to assess how the choice of genes 
#   (based on their information content) affects the inferred phylogenetic tree. 
# We select the top 90% (0.9 fraction) best scoring genes based on their information content metrics,
#   and infer a new phylogenetic tree using only these genes.
# Then we compare this tree to the reference tree by using the Robinson-Foulds metric, 
#    which quantifies the distance between two trees.
python3 gene-sensitivity-analysis.py \
    --input ${OUTPUT_DIR}/gene-metrics.tsv \
    --trimmed-alignments ${TRIMMED_ALIGNMENTS} \
    --output-dir ${OUTPUT_DIR}/gene-sensitivity-top-${TOP_FRACTION} \
    --top-fraction ${TOP_FRACTION} \
    --model ${MODEL} \
    --ref-tree ${OUTPUT_DIR}/SUPERMATRIX.aln.fasta.treefile \
    --threads ${THREADS}

conda deactivate
