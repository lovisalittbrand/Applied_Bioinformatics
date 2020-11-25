#!/bin/bash 
#===============================================================================
#
# FILE: bash_pipeline.sh
#
# DESCRIPTION: The following script performs and connects a series of bioinformatic analyses in a pipeline.
#			   The pipeline are performed in data located in folder /data
#			   The pipeline can run both constrained and unconstrained phylogenetic analyses with IQ-tree based on parameter input
#
# USAGE: There are two alternative ways of using this pipeline depending on if you have a concatenation or not.
# (1) $ sh bash_pipeline.sh COGS_FILE UNIQUE_TAXA MAPPING_ARCOGS
# (2) $ sh bash_pipeline.sh COGS_FILE UNIQUE_TAXA MAPPING_ARCOGS CONSTRAINT_TREE CONCAT_FILE PART_FILE
#
#	- COGS_FILE: input txt-file containing names of all arCOG-files that will be analyzed
#	- UNIQUE_TAXA: list containing names of all unique taxa used for mapping
#	- MAPPING_ARCOGS: mapping-file containing arCOG-names, gene names and taxa names 
#	- CONSTRAINT_TREE: tree-file containing the constrained topology
#	- CONCAT_FILE: concatenation-file of your aligned sequences
#	- PART_FILE: partition-file containing gene limits in your concatenated alignment
#
#===============================================================================

#-------------------------------------------------------------------------------
# DATA
#-------------------------------------------------------------------------------

cogs_file=$1
unique_taxa=$2
mapping_arcogs=$3
constraint_tree=$4
concat_file=$5
part_file=$6


cogs_vector=()
while read LINE  
do  
  cogs_vector+=($LINE)
done <$cogs_file


# --- Second run of pipeline: concatenation exist
if [ -f "$constraint_tree" ];
then
	#-------------------------------------------------------------------------------
	# PIPELINE STEP: 4		Creating phylogenetic tree
	#-------------------------------------------------------------------------------
	echo "\n#####"
	echo "Pipeline step 4: Creating phylogenetic tree based on constraint topology"
	echo "Data: Concatenation"
	echo "Software: IQ-Tree"
	echo "#####\n"
	
	iqtree -s $concat_file -spp $part_file -st "AA" -wsl -nt 1 -pre T2_ML -g $constraint_tree -m JTT+R3
	
	for f in T2_ML.*
	do
		mv $f "data/3_concatenation_iqtree/$f"
	done
	
# --- First run of pipeline: concatenation does not exist
else 
	#-------------------------------------------------------------------------------
	# PIPELINE STEP: 1-2	Alignment	Trimming
	#-------------------------------------------------------------------------------
	mkdir "data/1_alignment"
	mkdir "data/2_trimming"

	for i in ${cogs_vector[@]}
	do
		echo "\n#####"
		echo "Pipeline steps 1-3: Alignment; Trimming; Removing Line breaks"
		echo "Data: ${i}"
		echo "Software: MAFFT; Trimal-Gappyout"
		echo "#####\n"
	
		input="data/${i}.fasta"

		output_aln="data/1_alignment/${i}.mafft.fasta"
		mafft --auto $input > $output_aln


		output_trim="data/2_trimming/${i}.trim.mafft.fasta"
		trimal -gappyout -in $output_aln -out $output_trim

	done

	#-------------------------------------------------------------------------------
	# PIPELINE STEP: 3		Concatenation
	#-------------------------------------------------------------------------------

	echo "\n#####"
	echo "Pipeline step 3: Concatenation"
	echo "Data: ${cogs_file_list[@]}"
	echo "#####\n"

	mkdir "data/3_concatenation_iqtree"
	output_con="data/3_concatenation_iqtree/concatenation.fasta"

	python3 upp_concatArCOGalignments.py -t $unique_taxa -c $cogs_file -m $mapping_arcogs -s .trim.mafft.fasta -f data/2_trimming

	mv "concatenation.fasta" "data/3_concatenation_iqtree/concatenation.fasta"
	mv "partition.nexus" "data/3_concatenation_iqtree/partition.nexus"
	

	#-------------------------------------------------------------------------------
	# PIPELINE STEP: 4		Creating phylogenetic tree
	#-------------------------------------------------------------------------------

	echo "\n#####"
	echo "Pipeline step 4: Creating phylogenetic tree"
	echo "Data: Concatenation"
	echo "Software: IQ-Tree"
	echo "#####\n"

	partition_result="data/3_concatenation_iqtree/partition.nexus"
	concatenation_result="data/3_concatenation_iqtree/concatenation.fasta"


	iqtree -s $concatenation_result -spp $partition_result -st "AA" -wsl -nt 1 -pre T1_ML -m JTT+R3
	for f in T1_ML.*
	do
		mv $f "data/3_concatenation_iqtree/$f"
	done


	#-------------------------------------------------------------------------------
	# PIPELINE STEP: 5		Visualization frequency & entropy
	#-------------------------------------------------------------------------------

	echo "\n#####"
	echo "Pipeline step 5: Visualizing frequency and entropy"
	echo "Data: Concatenation"
	echo "#####\n"
	bin=1
	type="AA"

	python3 align_vis.py $concatenation_result $bin $type -pf $partition_result
	echo "Opening result in browser..."
fi
