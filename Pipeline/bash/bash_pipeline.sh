#!/bin/bash 
#===============================================================================
#
# FILE: bash_pipeline.sh
#
# DESCRIPTION: The following script performs and connects a series of bioinformatic analyses in a pipeline.
#			   The pipeline are performed in data located in folder /data
#
# USAGE: 
#  $ sh bash_pipeline.sh 
#
#===============================================================================

#-------------------------------------------------------------------------------
# DATA
#-------------------------------------------------------------------------------

#COG-files in FASTA-format that should be located in /data
cogs=(arCOG04070 arCOG04067)

#-------------------------------------------------------------------------------
# PIPELINE STEP: 1-2	Alignment	Trimming
#-------------------------------------------------------------------------------
mkdir "data/1_alignment"
mkdir "data/2_trimming"

for i in ${cogs[@]}
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

cogs_file_list=()
for y in ${cogs[@]} 
do 
	cogs_file_list+=(data/2_trimming/${y}.trim.mafft.fasta)
done

echo "\n#####"
echo "Pipeline step 3: Concatenation"
echo "Data: ${cogs_file_list[@]}"
echo "#####\n"

mkdir "data/3_concatenation"
output_con="data/3_concatenation/concatenation.fasta"

perl catfasta2phyml_mod.pl ${cogs_file_list[@]} -c -f -v > $output_con 


#-------------------------------------------------------------------------------
# PIPELINE STEP: 4		Removing line break
#-------------------------------------------------------------------------------

echo "\n#####"
echo "Pipeline step 4: Removing line break"
echo "Data: Concatenation"
echo "#####\n"

mkdir "data/4_linebreak"
output_lb="data/4_linebreak/concatenation.lb.fasta"
sh line_break.sh $output_con $output_lb

#-------------------------------------------------------------------------------
# PIPELINE STEP: 5		Creating phylogenetic tree
#-------------------------------------------------------------------------------

echo "\n#####"
echo "Pipeline step 5: Creating phylogenetic tree"
echo "Data: Concatenation"
echo "Software: IQ-Tree"
echo "#####\n"
mv "partition.nexus" "data/3_concatenation/partition.nexus"
partition_file="data/3_concatenation/partition.nexus"

#iqtree -s $output_con -spp $partition_file -st "AA" -wsl -nt 1

mkdir "data/5_iqtree"
mv "concatenation.fasta.log" "data/5_iqtree/concatenation.fasta.log"
mv "concatenation.fasta.bionj" "data/5_iqtree/concatenation.fasta.bionj"
mv "concatenation.fasta.ckp.gz" "data/5_iqtree/concatenation.fasta.ckp.gz"
mv "concatenation.fasta.mldist" "data/5_iqtree/concatenation.fasta.mldist"
mv "concatenation.fasta.model.gz" "data/5_iqtree/concatenation.fasta.model.gz"
mv "concatenation.fasta.treefile" "data/5_iqtree/concatenation.fasta.treefile"

#-------------------------------------------------------------------------------
# PIPELINE STEP: 6		Visualization frequency & entropy
#-------------------------------------------------------------------------------

echo "\n#####"
echo "Pipeline step 6: Visualizing frequency and entropy"
echo "Data: Concatenation"
echo "#####\n"
bin=10
type="AA"

python3 align_vis.py $output_lb $bin $type -pf $partition_file


#-------------------------------------------------------------------------------
# PIPELINE STEP: 7		Calculation & visualization phylogenetic signal 
#-------------------------------------------------------------------------------

echo "\n#####"
echo "Pipeline step 7: Calculating and visualizing phylogenetic signal"
echo "Data: Concatenation"
echo "#####\n"

slh="data/5_iqtree/concatenation.fasta.sitelh" 

#python3 SLS_diff.py $slh

#python3 GLS_diff.py $slh $partition_file

