#!/bin/bash 
#===============================================================================
#
# FILE: run_bash_p2.sh
#
# DESCRIPTION: The following script runs a part of the pipeline in "bash_pipeline.sh" and includes the following steps:
#	1. Phylogenetic tree reconstruction using tree constraint
#	It will thereafter combine files containing site likelihood from two trees and visualize phylogenetic signal.
#
# USAGE: 
#  $ sh run_bash_p2.sh 
#	
# REQUIREMENTS:
#	1. Txt-file containing a list of all arCOG-files (FASTA) that will be studied. 
#	2. File containing a list of all unique taxa names. 
#	3. Mapping-file containing a list of all arCOG, taxa and gene names. 
# 	4. Tree-file containing the constraint tree. 
#	5. Concatentation-file produced when running run_bash_p1.sh.
#	6. Partition-file produced when running run_bash_p1.sh.
#
# NOTE: Make sure paths to required input parameters to pipeline is correct.
#
#===============================================================================

#-------------------------------------------------------------------------------
# DATA
#-------------------------------------------------------------------------------

#Path to required parameters
cogs_file="data/cog_list.txt" 
unique_taxa="data/our_unique_taxa.list"
mapping_arcogs="data/mapping_arcogs.maps"
constraint_tree="data/T2_constraint.treefile"
concat_file="data/3_concatenation_iqtree/concatenation.fasta"
part_file="data/3_concatenation_iqtree/partition.nexus"

#-------------------------------------------------------------------------------
# RUN PIPELINE 
#-------------------------------------------------------------------------------

#Run to obtain site likelihoods from constraint tree 
sh bash_pipeline.sh --cogs $cogs_file --taxa $unique_taxa -m $mapping_arcogs --tree $constraint_tree --concat $concat_file -p $part_file

#-------------------------------------------------------------------------------
# MERGE SLH-FILES
#-------------------------------------------------------------------------------
#Create file containing site likelihood from unconstrained and constrained run.
awk '{sub(/Site_Lh/,"Tree1")}1' "data/3_concatenation_iqtree/T1_ML.sitelh" > "data/3_concatenation_iqtree/T1_ML_mod.sitelh"
awk '{sub(/Site_Lh/,"Tree2")}1' "data/3_concatenation_iqtree/T2_ML.sitelh" > "data/3_concatenation_iqtree/T2_ML_mod.sitelh"


touch "data/3_concatenation_iqtree/combined.sitelh"
sed -n 2p "data/3_concatenation_iqtree/T1_ML_mod.sitelh" >> "data/3_concatenation_iqtree/combined.sitelh"
sed -n 2p "data/3_concatenation_iqtree/T2_ML_mod.sitelh" >> "data/3_concatenation_iqtree/combined.sitelh"


#-------------------------------------------------------------------------------
# CALCULATION & VISUALIZATION: PHYLOGENETIC SIGNAL 
#-------------------------------------------------------------------------------
slh="data/3_concatenation_iqtree/combined.sitelh"
partition="data/3_concatenation_iqtree/partition.nexus"

#Calculate and visualize the phylogenetic signal for both individual sites and genes based on output from IQ-tree
#python3 SLS_diff.py $slh 
python3 GLS_diff.py $slh $partition -gs
