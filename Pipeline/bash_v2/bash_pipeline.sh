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
# (1) $ sh bash_pipeline.sh -c <COGS_FILE> --taxa <UNIQUE_TAXA> -m <MAPPING_ARCOGS>
# (2) $ sh bash_pipeline.sh -c <COGS_FILE> --taxa <UNIQUE_TAXA> -m <MAPPING_ARCOGS> --tree <CONSTRAINT_TREE> --concat <CONCAT_FILE> -p <PART_FILE>
#
#	- COGS_FILE: input txt-file containing names of all arCOG-files that will be analyzed
#	- UNIQUE_TAXA: list containing names of all unique taxa used for mapping
#	- MAPPING_ARCOGS: mapping-file containing arCOG-names, gene names and taxa names 
#	- CONSTRAINT_TREE: tree-file containing the constrained topology
#	- CONCAT_FILE: concatenation-file of your aligned sequences
#	- PART_FILE: partition-file containing gene limits in your concatenated alignment
#
# EXAMPLE: 
# (1) $ sh bash_pipeline.sh cog_list.txt unique_taxa.list mapping_arcogs.maps 
# (2) $ sh bash_pipeline.sh cog_list.txt unique_taxa.list mapping_arcogs.maps constraint.treefile concatenation.fasta partition.nexus 
#===============================================================================

#-------------------------------------------------------------------------------
# DATA
#-------------------------------------------------------------------------------

function usage() 
{
	echo ""
	echo "Usage: sh bash_pipeline [-c COGS_FILE] [-t TAXA_FILE] [-m MAPPING_FILE] [[-t TREE_CONSTRAINT] [--concat CONCATENATION] [-p PARTITION]]"
	echo ""
	echo "GENERAL OPTIONS:"
	echo "-h, --help	Printing help usages"
	echo ""
	echo "MANDATORY USAGE:"
	echo "-c, --cogs	Txt-file containing names of all arCOGS that will be studied"
	echo "--taxa		Txt-file containing all unique taxa names needed for mapping"
	echo "-m, --map	Mapping file containing all relationship"
	echo ""
	echo "OPTIONAL ADDITION OF PARAMETERS (ALL MUST BE INCLUDED):"
	echo "--tree		Tree-file containing constraint topology"
	echo "--concat	Concatenated alignment"
	echo "-p		Partition file (nexus)"
}



while [ "$1" != "" ]; do
    PARAM=`echo $1 #| awk -F= '{print $1}'`
    VALUE=`echo $2 #| awk -F= '{print $2}'`

	case $PARAM in 
		-c|--cogs)
		COGS_FILE=$VALUE
		;;
		--taxa)
		UNIQUE_TAXA=$VALUE
		;;
		-m|--map)
		MAPPING_COGS=$VALUE
		;;
		--tree)
		CONSTRAINT_TREE=$VALUE
		;;
		--concat)
		CONCAT_FILE=$VALUE
		;;
		-p)
		PART_FILE=$VALUE
		;;
		-h|--help)
		usage
		exit 1
		;;
		*)
		echo "ERROR: unknown parameter \"$PARAM\""		
		exit 1
		;;
	esac
	shift
	shift
done



cogs_vector=()
while read LINE  
do  
  cogs_vector+=($LINE)
done <$COGS_FILE


# --- Second run of pipeline: concatenation exist
if [ -f "$CONSTRAINT_TREE" ];
then
	#-------------------------------------------------------------------------------
	# PIPELINE STEP: 4		Creating phylogenetic tree
	#-------------------------------------------------------------------------------
	echo "\n#####"
	echo "Pipeline step 4: Creating phylogenetic tree based on constraint topology"
	echo "Data: Concatenation"
	echo "Software: IQ-Tree"
	echo "#####\n"
	
	iqtree -s $CONCAT_FILE -spp $PART_FILE -st "AA" -wsl -nt 1 -pre T2_ML -g $CONSTRAINT_TREE -m JTT+R3
	
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

	python3 upp_concatArCOGalignments.py -t $UNIQUE_TAXA -c $COGS_FILE -m $MAPPING_COGS -s .trim.mafft.fasta -f data/2_trimming

	mv "concatenation.fasta" "data/3_concatenation_iqtree/concatenation.fasta"
	mv "partition.nexus" "data/3_concatenation_iqtree/partition.nexus"
	

	#-------------------------------------------------------------------------------
	# PIPELINE STEP: 4		Creating phylogenetic tree
	#-------------------------------------------------------------------------------

	echo "\n#####"
	echo "Pipeline step 4: Creating phylogenetic tree"
	echo "Data: Concatenation"
	echo "Software: Iqtree"
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
	bin=10
	type="AA"

	python3 align_vis.py $concatenation_result $bin $type -pf $partition_result
	echo "Opening result in browser..."
fi
