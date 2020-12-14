This folder contains the scripts of the project, used in the pipeline to obtain the log-likelihoods and to make the visualization plots.
* aa_color.txt - contains the coloring scheme of the amino acids in the frequency plot. It can be modified according to the user's preference.
* align_vis.py - used to visualise the base frequency from a MSA. It creates plots for amino acid and nucleotide frequency per site, amino acid frequency per taxa, entropy and a gene localization bar.
* catfasta2phyml.pl - used to concatenate aligned sequences that will be then used to reconstruct a phylogeny
* catfasta2phyml_mod.pl - same as above but this version was modified to print the partition file
* dna_color.txt - contains the coloring scheme of the nucleotides in the frequency plot. It can be modified according to the user's preference.
* signal_vis2trees.py - visualize phylogenetic signal per site and per gene for two conflicting topologies.
* signal_vis3trees.py - visualize phylogenetic signal per site and per gene for three conflicting topologies
* site_remover.py - removes genes in order to investigate the effect of their absence on gene signal, support and tree topology.
