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
# 	1. Tree-file containing the constraint tree. 
#	2. Concatentation-file produced when running run_bash_p1.sh.
#
# OPTIONAL ADDITIONS:
#	1. Partition-file produced when running run_bash_p1.sh, matching the input concatenation, in case wanting to do partitioned phylogenetic reconstruction. 
#
# NOTE: Make sure paths to required input parameters to pipeline is correct.
#
#===============================================================================

#-------------------------------------------------------------------------------
# DATA
#-------------------------------------------------------------------------------

#Path to required parameters
constraint_tree="data/T2_constraint.treefile"
concat_file="data/3_concatenation/concatenation.fasta"
part_file="data/3_concatenation/partition.nexus"

#-------------------------------------------------------------------------------
# RUN PIPELINE 
#-------------------------------------------------------------------------------

#Run to obtain site likelihoods from constraint tree 
sh bash_pipeline.sh --tree $constraint_tree --concat $concat_file -p $part_file

#-------------------------------------------------------------------------------
# MERGE SLH-FILES
#-------------------------------------------------------------------------------
#Create file containing site likelihood from unconstrained and constrained run.
awk '{sub(/Site_Lh/,"Tree1")}1' "data/4_iqtree/T1_ML.sitelh" > "data/4_iqtree/T1_ML_mod.sitelh"
awk '{sub(/Site_Lh/,"Tree2")}1' "data/4_iqtree/T2_ML.sitelh" > "data/4_iqtree/T2_ML_mod.sitelh"


touch "data/4_iqtree/combined.sitelh"
echo "This is a file containing site likelihoods from 2 topologies" >> "data/4_iqtree/combined.sitelh"
sed -n 2p "data/4_iqtree/T1_ML_mod.sitelh" >> "data/4_iqtree/combined.sitelh"
sed -n 2p "data/4_iqtree/T2_ML_mod.sitelh" >> "data/4_iqtree/combined.sitelh"


#-------------------------------------------------------------------------------
# CALCULATION & VISUALIZATION: PHYLOGENETIC SIGNAL 
#-------------------------------------------------------------------------------
slh="data/4_iqtree/combined.sitelh"
partition="data/3_concatenation/partition.nexus"

#Calculate and visualize the phylogenetic signal for both individual sites and genes based on output from IQ-tree 
python3 signal_vis2trees.py $slh $partition -gs -gt -ss
