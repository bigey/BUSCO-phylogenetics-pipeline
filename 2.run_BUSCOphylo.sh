#!/bin/bash
set -euo pipefail

BUSCO_RESULTS="Busco-results"
OUT_PHYLO_DIR="Phylogenomics"
THREADS=64

conda activate busco-phylo

# It is sometimes required to unset the MAFFT_BINARIES environment variable to allow BUSCOphylo.py to find the correct MAFFT executable
unset MAFFT_BINARIES

# Run the BUSCOphylo.py script to generate the supermatrix, supermatrix phylogeny and individual supertrees
#   --supermatrix option will generate the concatenated alignment phylogeny
#   --supertree option will generate individual gene trees and then concatenate them
#   --concordance option is used to calculate gene concordance factors for the supermatrix phylogeny

python3 BUSCOphylo.py --directory ${BUSCO_RESULTS} --output ${OUT_PHYLO_DIR} --supermatrix --supertree --concordance --threads ${THREADS}

conda deactivate