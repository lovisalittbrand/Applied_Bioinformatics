#!/bin/bash
#Script executing Snakemake pipeline


## Mandatory parameters:
# DATA = Specifying name of input file: DATA.fasta
# TYPE = DNA or protein sequence as input_file

## Optional parameters:
# BB = number of bootstrap replicates
# MODEL = substitution model if creating a phylogeny
# BIN = Bin size when visualizing amino acid frequency


## Visualizing DAG-graphs
DATA="zika" snakemake --cores 1 --dag phylo | dot -Tsvg > dag_phylo.svg
#DATA="zika" snakemake --cores 1 --dag plot_freq | dot -Tsvg > dag_plot.svg

## Running pipeline with different options
#DATA="zika" BB=1000 TYPE="DNA" MODEL="GTR+F+I+G4" snakemake -s Snakefile --cores 1 phylo
#Type="DNA" BIN=100 snakemake --cores 1 plot_freq
#snakemake -s data/Snakefile --cores 1 phylogenetic_signal

#TYPE="AA" snakemake -s Snakefile --cores 1 kasia_data/arCOG00779.lb.trim.mafft.fasta kasia_data/arCOG00987.lb.trim.mafft.fasta
TYPE="AA" snakemake -s Snakefile --cores 1 kasia_data/concat.lb.trim.mafft.fasta
