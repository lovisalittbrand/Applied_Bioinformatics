#!/bin/bash 
#===============================================================================
#
# FILE: run_bash_p1.sh
#
# DESCRIPTION: The following script runs the pipeline in "bash_pipeline.sh" and includes the following steps:
#	1. Align 
#	2. Trim 
#	3. Concatenation
#	4. Phylogenetic tree reconstruction
#	5. Frequency & Entropy visualization of concatenation
#
# USAGE: 
#  $ sh run_bash_p1.sh 
#
# REQUIREMENTS:
#	1. Txt-file containing a list of all arCOG-files (FASTA) that will be studied. 
#	2. File containing a list of all unique taxa names. 
#	3. Mapping-file containing a list of all arCOG, taxa and gene names. 
#
# NOTE: Make sure paths to required input parameters to pipeline is correct.
#	
#===============================================================================

#-------------------------------------------------------------------------------
# DATA
#-------------------------------------------------o------------------------------

#Path to required parameters
cogs_file="data/25AERibosomalProteins.lines" 
unique_taxa="data/AERibosomalProteins.taxa"
mapping_arcogs="data/AERibosomalProteins.maps"


#-------------------------------------------------------------------------------
# RUN PIPELINE
#-------------------------------------------------------------------------------

#Run to obtain site likelihoods from unconstraint tree 
sh bash_pipeline.sh --cogs $cogs_file --taxa $unique_taxa -m $mapping_arcogs