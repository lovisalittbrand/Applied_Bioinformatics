#!/bin/bash 
#===============================================================================
#
# FILE: run_bash.sh
#
# DESCRIPTION: The following script runs the pipeline in "bash_pipeline.sh" twice with different tree constraints.
#			   The script thereafter visualize the phylogenetic signal per site and gene
#
# USAGE: 
#  $ sh run_bash.sh 
#	
#===============================================================================

#-------------------------------------------------------------------------------
# DATA
#-------------------------------------------------------------------------------

#Path to tree constraint file and txt-file listing all arCOGs
tree_constraint="data/partition.nexus.treefile"
cogs_file="data/cog_list.txt" 

#-------------------------------------------------------------------------------
# RUN PIPELINE
#-------------------------------------------------------------------------------

#Run to obtain unconstraint tree 
sh bash_pipeline.sh $cogs_file

#Run to obtain constraint tree
sh bash_pipeline.sh $cogs_file $tree_constraint

#-------------------------------------------------------------------------------
# COMPUTE SITE LIKELIHOOD FOR UNCONSTRAINED VS. CONSTRAINED TREE
#-------------------------------------------------------------------------------

#Create trees-file with both constrained and unconstrained tree
touch "data/4_linebreak_iqtree/trees.tree"
echo "$(cat data/4_linebreak_iqtree/T1_ML)" >> "data/4_linebreak_iqtree/trees.tree" #CHANGE PATHS!!!!
echo "$(cat data/4_linebreak_iqtree/T2_ML)" >> "data/4_linebreak_iqtree/trees.tree"

#Compute site likelihood using both trees
iqtree -wsl -nt 1 -st AA -z "data/4_linebreak_iqtree/trees.tre" -s "data/4_linebreak_iqtree/concatenation.lb.fasta" -spp "data/4_linebreak_iqtree/partition.nexus" -pre sum_site_lh

#-------------------------------------------------------------------------------
# CALCULATION & VISUALIZATION: PHYLOGENETIC SIGNAL 
#-------------------------------------------------------------------------------
slh="data/4_lineabreak_iqtree/sum_site_lh.sitelh"
partition="data/4_linebreak_iqtree/partition.nexus"

#Calculate and visualize the phylogenetic signal for both individual sites and genes based on output from IQ-tree
python3 SLS_diff.py $slh 
python3 GLS_diff.py $slh $partition