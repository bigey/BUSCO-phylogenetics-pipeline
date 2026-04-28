# BUSCO Phylogenomics pipeline

## Presentation

The following scripts are used to construct species phylogenies using BUSCO single copy proteins. It works directly from BUSCO outputs and can be used for supermatrix or supertree/coalescent methods. The program will automatically identify single-copy BUSCO proteins, generate alignments using `MAFFT` and `ClipKIT`. Then it either concatenates them into a supermatrix fasta fileto infer the species-tree phylogeny using `IQ-TREE` or generate individual trees for a supertree approach. The program can also perform gene and sequence concordance factors (gCF and sCF) analysis if both supermatrix and supertree methods are selected (`--concordance`). The resulting supermatrix species tree (in newick format) is labeled with gene and sequence concordance factors (`gCF` and `sCF`). This can provide insights into the level of gene tree discordance and the robustness of the inferred species tree. 

The pipeline is designed to perform sensitivity analysis by evaluating the impact of the number of genes included in the analysis on the resulting species tree topology. It computes several metrics for each gene tree (alignment length, average bipartition support, relative composition variability, median long branch score, treeness, saturation, and treeness/RCV ratio) and subsets the genes based on the specified metric and fraction (e.g., top 25% of genes based on alignment length). Then it infers new trees using the subseted genes and compares the resulting trees to a reference tree (in this case, the supermatrix tree) using the Robinson-Foulds distance metric. This allows to evaluate how the number of genes included in the analysis impacts the resulting species tree topology.

## Version

+ v1.0 (2026-03-08): Initial release, including supermatrix and supertree methods, concordance factor computation.
+ v2.0 (2026-04-28): Added gene sensitivity analysis, including gene metrics computation and tree comparison using Robinson-Foulds distance.

## History

Script `busco-phylo.py` is derived from the original pipeline available at [GitHub](https://github.com/jamiemcg/BUSCO_phylogenomics.git), which was developed by [Jamie McGowan](https://jamiemcgowan.ie/) and is licensed under the MIT License. It has been modified to include additional features and improvements like concordance factor analysis. 

The sensitivity analysis part of the pipeline (`compute-gene-metrics.py` and `gene-sensitivity-analysis.py`) is inspired by the sensitivity analysis tutorial of [Jacob L. Steenwyk](https://jlsteenwyk.com). Specifically, the gene sensitivity analysis is based on the tutorial available at [Sensitivity analysis: detecting incongruence at different scales](https://jlsteenwyk.com/tutorials/ub_sensitivity_analysis.html).

## Requirements

The following softwares and packages should be installed to run the pipeline:

* [Python3](https://www.python.org/)
* [BioPython](https://biopython.org/)
* [pandas](https://pandas.pydata.org/)
* [matplotlib](https://matplotlib.org/)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [ClipKIT](https://jlsteenwyk.com/ClipKIT/)
* [PhyKIT](https://jlsteenwyk.com/PhyKIT/)
* [IQ-TREE](https://iqtree.github.io/)

The simple way is to create a conda environment where all the dependencies are installed. A `environment.yml` file is provided for this purpose. You can create a dedicated conda environment with the following command:

```shell
conda env create -f environment.yml
```

## Principle of the pipeline

For each species, BUSCO creates a directory containing all the results needed to perform the analysis. Inside this directory, there are several sub-directories, but the one we are interested in is `run_[lineage]_odb[XX]`, which contains the single-copy BUSCO protein sequences in fasta format. In order to centralize the all the BUSCO files for multiple species analysis, create a new directory (`BUSCO_RESULTS`) and link all the `run_[lineage]_odb[XX]` sub-directories to this directory. Link names should start by `run_`, then folowed by the name of the species `species1`.

Exemple of the directory structure:

```shell
$ tree BUSCO_RESULTS/
BUSCO_RESULTS/
├── run_species1 -> /path/to/BUSCO_dir_species1/run_[lineage]_odb[XX]
├── run_species2 -> /path/to/BUSCO_dir_species2/run_[lineage]_odb[XX]
├── run_species3 -> /path/to/BUSCO_dir_species3/run_[lineage]_odb[XX]
├── run_species4 -> /path/to/BUSCO_dir_species4/run_[lineage]_odb[XX]
├── run_species5 -> /path/to/BUSCO_dir_species5/run_[lineage]_odb[XX]
├── run_species6 -> /path/to/BUSCO_dir_species6/run_[lineage]_odb[XX]
```

Then you can use this directory as the input for `busco-phylo.py` script with the `--directory` option.

The utility script `1.run_busco_genome_mode.sh` can be used to automate running BUSCO for a set of genomes present in a given directory (in fasta format, *.fasta). See below for more details on how to use this script.

### Supermatrix method

For exemple, to perform phylogenetic inference using the supermatrix method, use the following command:

```shell
python3 busco-phylo.py --directory BUSCO_RESULTS --output OUTPUT_DIRECTORY --supermatrix --threads 8
```

This will generate a concatenated alignment of all single-copy BUSCO proteins (using `MAFFT` and `ClipKIT`) and perform phylogenetic inference using `IQ-TREE`. The resulting alignment and tree will be saved in the output directory (`OUTPUT_DIRECTORY`).

You can use the `--stop-early` option which stops the pipeline just after generating the concatenated alignment (supermatrix). Then you can take a look at the concatenated alignment and manually choose parameters for phylogenetic inference.

### Supertree method

To perform a supertree analysis, use the following command:

```shell
python3 busco-phylo.py --directory BUSCO_RESULTS --output OUTPUT_DIRECTORY --supertree --threads 8
```

This command performs alignment and phylogenetic inference for each single-copy BUSCO protein family separately and then combines the individual trees into a  concatenated species trees file `ALL.trees`. All the data will be saved in the output directory (`OUTPUT_DIRECTORY`).

A species tree could be inferred using a coalescent-based approach that can account for gene tree discordance due to incomplete lineage sorting.

### Supermatrix and supertree methods plus concordance analysis

To run both supermatrix and supertree methods and perform concordance factor analysis (gCF and sCF), use the following command:

```shell
python3 busco-phylo.py --directory BUSCO_RESULTS --output OUTPUT_DIRECTORY --supermatrix --supertree --concordance --threads 8
```

The resulting supermatrix species tree (in newick format) is labeled with gene and sequence concordance factors (`gCF` and `sCF`). This can provide insights into the level of gene tree discordance and the robustness of the inferred species tree.

At the end of the run, you will find the resulting alignments and trees in the `OUTPUT_DIRECTORY` directory:

* `SUPERMATRIX.aln.fasta`: the concatenated supermatrix alignment,
* `ALL.trees`: the concatenated family protein trees,
* `SUPERMATRIX.treefile`: the supermatrix tree with bootstrap support values (in newick format),
* `gcf.cf.tree`: the supermatrix tree labeled with gene concordance factor (gCF),
* `sCF.cf.tree`: the supermatrix tree labeled with sequence concordance factor (sCF).

Branch labels in the `SUPERMATRIX.treefile.cf.tree` file are formatted as follows:

`(species1:branch_length,species2:branch_length)bootstrap/gCF/sCF:branch_length,`

Where `bootstrap` is the bootstrap support value for the branch, `gCF` is the gene concordance factor, and `sCF` is the sequence concordance factor (percentage). Values are separated by slashes (`/`) and the branch length (float) is indicated after the colon (`:`).

### `busco-phylo.py` required parameters

* `--directory`: input directory containing BUSCO results (each link should be named `run_[species]` and point to the corresponding BUSCO `run_[lineage]_odb[XX]` sub-directory)
* `--output`: phylogenetic output directory
* `--supermatrix`: choose to run supermatrix method
* `--supertree`: choose to run supertree method

### `busco-phylo.py` optional parameters

* `--percent-single-copy`: fraction of the BUSCO single-copy genes species will be included in the analysis (default = 1.0).
* `--model`: protein evolution model to use for IQ-TREE phylogeny inference (see documentation, default to LG+R4+F)".
* `--concordance`: run gene concordance factor analysis (automatic turns on `--supermatrix` and `--supertree` methods).
* `--stop-early`: stop pipeline after generating concatenated alignment, only relevant with `--supermatrix` method, incompatible with `--supertree` and `--concordance`.
* `--threads`: number of threads to use (default = 8)

## Running BUSCO for multiple species (step 1)

An utility script is provided to automate the runs of BUSCO for a set of genomes present in a given directory (in fasta format, *.fasta). It generates the required input directory structure for the phylogenomics pipeline.

Adapt the parameters in `1.run_busco_genome_mode.sh`:

```shell
GENOME_DIR="Genomes"           # Directory containing the sequences (fasta), should be named `species1.fasta`, `species2.fasta`, etc.
BUSCO_RESULTS="Busco-results"  # Directory where BUSCO results will be centralized, will contain symbolic links to the per-species BUSCO result directories
LINEAGE="saccharomycetes"      # Set the appropriate lineage for your species, see BUSCO documentation for available lineages, (or command line: busco --list-datasets)
PREDICTOR="metaeuk"            # Select an appropriate gene predictor for BUSCO genome mode (options: miniprot, metaeuk, augustus)
THREADS=16                     # Number of threads to use for BUSCO runs
```

Then run the script:

```shell
bash -i 1.run_busco_genome_mode.sh
```

The script will run BUSCO for each genome in `GENOME_DIR` and create symbolic links to the BUSCO result directories in `BUSCO_RESULTS`. You can then use `--directory BUSCO_RESULTS` to run the the phylogenomics script `busco-phylo.py`.

## Running the phylogenomics pipeline (step 2)

The second helper script `2.run_busco_phylogenomics.sh` is a wrapper to run the phylogenomics pipeline with the desired parameters.

Adapt the parameters in the script:

```shell
BUSCO_RESULTS="Busco-results" # Directory containing the centralized BUSCO results (with links to the per-species BUSCO result directories)
OUT_PHYLO_DIR="Phylogenetics" # Directory where phylogenetic results will be saved (see above for the output files generated by the pipeline)
MODEL="LG+R4+F"               # Protein evolution model to use for IQ-TREE phylogeny inference (see IQ-TREE documentation for available models)
THREADS=16                    # Number of threads you want to allocate for the phylogenomic pipeline
```

Then run the script:

```shell
bash -i 2.run_busco-phylogenomics.sh
```

## Performing gene sensitivity analysis (step 3)

The third step involves performing a gene sensitivity analysis to evaluate the impact of the number of genes included in the analysis on the resulting species tree topology. The script `3.run_gene_sensitivity_analysis.sh` is a wrapper to run the gene sensitivity analysis with the desired parameters.

Adapt the parameters in the script:

```shell
TRIMMED_ALIGNMENTS="Phylogenetics/trimmed-alignments" # Directory containing the trimmed alignments for each gene (output of the phylogenomics pipeline, step 2)
TREES_DIR="Phylogenetics/trees"                       # Directory containing the gene trees for each gene (output of the phylogenomics pipeline, step 2)
OUTPUT_DIR="Gene-sensitivity-analysis"                # Directory where gene sensitivity analysis results will be saved
TOP_FRACTION=0.25                                     # Set the fraction of genes to include in the analysis (e.g., 0.25 for top 25% of best scoring genes)
MODEL="LG+R4+F"                                       # Set the protein evolution model to use for IQ-TREE phylogeny inference (see documentation for available models)
THREADS=64                                            # Set the number of threads you want to allocate for the phylogenomic pipeline
```

Then run the script:

```shell
bash -i 3.run_gene_sensitivity_analysis.sh
```
The script will first compute several metrics for each gene tree (alignment length, average bipartition support, relative composition variability, median long branch score, treeness, saturation, and treeness/RCV ratio). Then it will subset the genes based on the specified metric and fraction (e.g., top 25% of genes based on alignment length), infer new trees using the subseted genes, and compare the resulting trees to the reference tree (in this case, the supermatrix tree) using the Robinson-Foulds distance metric. The results will be saved in the specified output directory for further analysis.

### First, compute gene metrics

First, the script computes several metrics for each gene tree. This script will read the trimmed alignments and their corresponding gene trees, and compute several metrics. Specifically, it will calculate:

+ Alignment length (`aln_len`): longer alignments are associated with greater phylogenetic signal
+ Average bipartition support (`abs`): single-locus phylogenies that are better supported are thought to have more robust phylogenetic signal
+ Relative composition variability (`rcv`): a measure of compositional bias wherein lower values are indicative of less bias
+ Median long branch score (`lbs`): lower scores are more desirable because they are indicative of phylogenetic trees less susceptible to long branch attraction artifacts
+ Treeness (`treeness`): a measure of the amount of phylogenetic signal in a gene tree; higher values are more desirable because they indicate that a greater proportion of the total tree length is attributable to internal branches, which are more likely to reflect true evolutionary relationships rather than noise or artifacts.
+ Saturation (`saturation`): a measure of the accuracy between distances measured using alignments compared to phylogenetic trees; data with no saturation will have a value of 1
+ Treeness/RCV ratio (`treeness_over_rcv`): Higher treeness/RCV values are thought to be desirable because they have a higher signal-to-noise ratio and are less susceptible to compositional biases.

```shell
python3 compute-gene-metrics.py \
  --trimmed-alignments ${TRIMMED_ALIGNMENTS} \
  --trees ${TREES_DIR} \
  --output ${OUTPUT_DIR}/gene-metrics.tsv \
  --threads ${THREADS}
````

The results are saved in a TSV file for further analysis (`${OUTPUT_DIR}/gene-metrics.tsv`).

This file contains three columns: `gene`, `metric`, and `value`. The `gene` column contains the name of the gene, the `metric` column contains the name of the metric (e.g., `aln_len`, `abs`, etc.), and the `value` column contains the corresponding value for that metric.

The file is used to plot the distribution of each metric across all genes (script `plot-gene-metrics.py`), and to perform gene sensitivity analysis by subsetting genes based on their metric values (e.g., top 25% of genes based on alignment length, etc.).

### Second, perform gene sensitivity analysis

Second, the script performs a gene sensitivity analysis by inferring species trees using subsets of genes. This allows to evaluate how the number of genes included in the analysis impacts the resulting species tree topology. For each metric, the genes are subseted based on the `--top-fraction` parameter (e.g., top 0.25, top 0.5, etc.). The script will read the trimmed alignments for the subseted genes, and infer new trees. The resulting trees are then compared to the reference tree (in this case, the supermatrix tree, `--ref-tree`) using the Robinson-Foulds distance metric. This allows to evaluate how the number of genes included in the analysis impacts the resulting species tree topology.

```shell
python3 gene-sensitivity-analysis.py \
  --input ${OUTPUT_DIR}/gene-metrics.tsv \
  --trimmed-alignments ${TRIMMED_ALIGNMENTS} \
  --output-dir ${OUTPUT_DIR}/gene-sensitivity-top-${TOP_FRACTION} \
  --top-fraction ${TOP_FRACTION} \
  --model ${MODEL} \
  --ref-tree ${OUTPUT_DIR}/SUPERMATRIX.aln.fasta.treefile \
  --threads ${THREADS}
```

All resulting files are saved in a directory (`--output-dir`) for further analysis. The main output files include:

+ `rf_distances.tsv`: a TSV file containing the Robinson-Foulds distances between the reference tree and the trees inferred using subsets of genes.
+ `<metric>.treefile`: the resulting trees inferred using subsets of genes based on the specified metric (e.g., `aln_len.treefile`, `abs.treefile`, etc.).

## Citation

If you use this code in your research, please cite:

BibTeX

```bibtex
@misc{bigey2026,
  author       = {Bigey, Frédéric},
  title        = {BUSCO Phylogenomics pipeline},
  year         = {2026},
  howpublished = {\url{https://github.com/bigey/BUSCO-phylogenetics-pipeline}},
  note         = {accessed 2026-XX-XX}
}
```

Biblatex

```biblatex
@software{bigey2026,
  author       = {Bigey, Frédéric},
  title        = {BUSCO Phylogenomics pipeline},
  year         = {2026},
  url          = {https://github.com/bigey/BUSCO-phylogenetics-pipeline},
  note         = {accessed 2026-XX-XX}
}
```
