# BUSCO Phylogenomics pipeline

## Presentation

Utility scripts to construct species phylogenies using BUSCO single copy genes. Works directly from BUSCO outputs and can be used for supermatrix or supertree/coalescent methods. The program will automatically identify single-copy BUSCO proteins, generate alignments (`MAFFT`) and trim (`CLipKit`) them. Then it either concatenates them into a supermatrix to infer the species-tree phylogeny using `IQ-TREE` or generate individual gene trees for a supertree approach. The program can also perform gene and sequence concordance factor (gCF and sCF) analysis if both supermatrix and supertree methods are selected (`--concordance`).

## History

This is a fork of the original pipeline available at [GitHub](https://github.com/jamiemcg/BUSCO_phylogenomics.git), which was developed by [Jamie McGowan](https://github.com/jamiemcg/BUSCO_phylogenomics.git) and is licensed under the MIT License.

It has been modified to include additional features and improvements like concordance factor analysis. Utility scripts have been added to automate the runs of BUSCO for multiple species and to run the phylogenomics pipeline with the desired parameters.

## Requirements

The following softwares and packages should be installed to run the pipeline:

* [Python3](https://www.python.org/)
* [BioPython](https://biopython.org/)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [ClipKit](https://jlsteenwyk.com/ClipKIT/)
* [IQ-TREE](https://iqtree.github.io/)

The simple way is to create a conda environment where all the dependencies are installed. A `environment.yml` file is provided for this purpose. You can create the conda environment with the following command:

```shell
conda env create -f environment.yml
```

## Usage

For each species, BUSCO creates a directory containing all the results needed to perform the analysisa. Inside this directory, there are several sub-directories, but the one we are interested in is `run_[lineage]_odb[XX]`, which contains the single-copy BUSCO protein sequences in fasta format. In order to centralize the BUSCO results for multiple species, create a new directory (`BUSCO_RESULTS`) and link all the `run_[lineage]_odb[XX]` sub-directories to this directory. The link name should start by `run_` then folowed by the name of the species.

This directory should have the following structure and naming convention:

* run_[species1] -> link to BUSCO `run_[lineage]_odb[XX]` sub-directory for species1
* run_[species2] -> link to BUSCO `run_[lineage]_odb[XX]` sub-directory for species2
* run_[species3] -> link to BUSCO `run_[lineage]_odb[XX]` sub-directory for species3
* run_[species4] -> link to BUSCO `run_[lineage]_odb[XX]` sub-directory for species4
* run_[species5] -> link to BUSCO `run_[lineage]_odb[XX]` sub-directory for species5
* run_[species6] -> link to BUSCO `run_[lineage]_odb[XX]` sub-directory for species6
* etc.

Then you can use this directory as the input for the phylogenomics pipeline.

The utility script `1.run_busco_genome_mode.sh` can be used to automate running BUSCO for a list of species, given a directory containing all the genome sequences (fasta format, *.fasta). It generates the required directory structure for the phylogenomics pipeline.

### Supermatrix method

For exemple, to perform phylogenetic inference using the supermatrix method, use the following command:

```shell
python3 BUSCOphylo.py --directory BUSCO_RESULTS --output OUTPUT_DIRECTORY --supermatrix --threads 8
```

This will generate a concatenated alignment of all single-copy BUSCO proteins (using `MAFFT` and `ClipKit`) and perform phylogenetic inference using `IQ-TREE`. The resulting alignment and tree will be saved in the output directory (`OUTPUT_DIRECTORY`).

You can use the `--stop-early` option which stops the pipeline just after generating the concatenated alignment (supermatrix). Then you can take a look at the concatenated alignment and manually choose parameters for phylogenetic inference.

### Supertree method

To perform a supertree analysis, use the following command:

```shell
python3 BUSCOphylo.py --directory BUSCO_RESULTS --output OUTPUT_DIRECTORY --supertree --threads 8
```

This command performs alignment and phylogenetic inference for each single-copy BUSCO protein family separately and then combines the individual trees into a  concatenated species trees file `ALL.trees`. All the data will be saved in the output directory (`OUTPUT_DIRECTORY`).

A species tree could be inferred using a coalescent-based approach that can account for gene tree discordance due to incomplete lineage sorting.

### Supermatrix and supertree methods plus concordance analysis

To run both supermatrix and supertree methods and perform concordance factor analysis (gCF and sCF), use the following command:

```shell
python3 BUSCOphylo.py --directory BUSCO_RESULTS --output OUTPUT_DIRECTORY --supermatrix --supertree --concordance --threads 8
```

The resulting supermatrix species tree is labeled with gene concordance factors (gCF and sCF), which represent the proportion of gene trees that support each branch in the species tree. This can provide insights into the level of gene tree discordance and the robustness of the inferred species tree.

### `BUSCOphylo.py` required parameters

* `--directory`: input directory containing BUSCO results (each link should be named `run_[species]` and point to the corresponding BUSCO `run_[lineage]_odb[XX]` sub-directory)
* `--output`: phylogenetic output directory
* `--supermatrix`: choose to run supermatrix method
* `--supertree`: choose to run supertree method

### `BUSCOphylo.py` optional parameters

* `--percent-single-copy`: fraction of the BUSCO single-copy genes species will be included in the analysis (default = 1.0).
* `--model`: protein evolution model to use for IQ-TREE phylogeny inference (see IQ-TREE documentation, default to LG method)".
* `--concordance`: run gene concordance factor analysis (automatic turns on `--supermatrix` and `--supertree` methods).
* `--stop-early`: stop pipeline after generating concatenated alignment, only relevant with `--supermatrix` method, incompatible with `--supertree` and `--concordance`.
* `--threads`: number of threads to use (default = 8)

## Running BUSCO for multiple species

An utility script is provided to automate the runs of BUSCO for a list of species, given a directory containing all the genome sequences (fasta format). It generates the required input directory structure for the phylogenomics pipeline.

Adapt the parameters in `1.run_busco_genome_mode.sh`:

```shell
GENOME_DIR="Genomes"           # Directory containing genome sequences in fasta format, should be named `species1.fasta`, `species2.fasta`, etc.
BUSCO_RESULTS="Busco-results"  # Directory where BUSCO results will be centralized, will contain symbolic links to the per-species BUSCO result directories
LINEAGE="saccharomycetes"      # Set the appropriate lineage for your species, see BUSCO documentation for available lineages, (or command line: busco --list-datasets)
PREDICTOR="metaeuk"            # Select an appropriate gene predictor for BUSCO genome mode, options are: miniprot, metaeuk, augustus.
THREADS=16                     # Number of threads to use for BUSCO runs
```

Then run the script:

```shell
bash 1.run_busco_genome_mode.sh
```

The script will run BUSCO for each genome in `GENOME_DIR` and create symbolic links to the BUSCO result directories in `BUSCO_RESULTS`. You can then use `BUSCO_RESULTS` as the input directory for the phylogenomics pipeline.

## Running the phylogenomics pipeline

The second script `2.run_phylo_pipeline.sh` is a wrapper to run the phylogenomics pipeline with the desired parameters. Adapt the parameters in the script:

```shell
BUSCO_RESULTS="Busco-results" # Directory containing the centralized BUSCO results, should be the same as the output from the previous script
OUT_PHYLO_DIR="Phylogenetics" # Directory where phylogenetic results will be saved
```

Then run the script:

```shell
bash 2.run_phylo_pipeline.sh
```

## Citation

If you use this document or script in your research, please cite:

BibTeX

```bibtex
@misc{bigey2026,
  author       = {Bigey, Frédéric},
  title        = {BUSCO Phylogenomics pipeline},
  year         = {2026},
  howpublished = {\url{https://github.com/bigey/busco-phylo}},
  note         = {accessed 2026-03-08}
}
```

Biblatex

```biblatex
@software{bigey2026,
  author       = {Bigey, Frédéric},
  title        = {BUSCO Phylogenomics pipeline},
  year         = {2026},
  version      = {v1.0.0},
  url          = {https://github.com/bigey/busco-phylo},
  note         = {accessed 2026-03-08}
}
```
