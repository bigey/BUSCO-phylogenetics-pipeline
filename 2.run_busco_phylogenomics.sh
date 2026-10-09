#!/bin/bash
set -euo pipefail

BUSCO_RESULTS="BUSCO-results"
OUT_PHYLO_DIR="Phylogenomics"
MODEL="LG+R4+F"
OUTGROUP="DebSin_CBS10405"
THREADS=64

conda activate env-busco-phylo

# It is sometimes required to unset the MAFFT_BINARIES environment variable to allow BUSCOphylo.py to find the correct MAFFT executable
unset MAFFT_BINARIES

# Run the BUSCOphylo.py script to generate the supermatrix, supermatrix phylogeny and individual supertrees
#   --supermatrix option will generate the concatenated alignment phylogeny
#   --supertree option will generate individual gene trees and then concatenate them
#   --concordance option is used to calculate gene concordance factors for the supermatrix phylogeny
#   --outgroup is used to "root" the tree (optional)

python3 busco-phylo.py \
    --busco-dir ${BUSCO_RESULTS} \
    --output-dir ${OUT_PHYLO_DIR} \
    --model ${MODEL} \
    --outgroup ${OUTGROUP} \
    --threads ${THREADS} \
    --supermatrix \
    --supertree \
    --concordance

conda deactivate
