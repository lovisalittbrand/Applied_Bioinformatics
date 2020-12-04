#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J LovisaLittbrand_Applied_Bioinformatics
#SBATCH --mail-type=ALL
#SBATCH --mail-user lollo.littbrand@hotmail.com

#Load appropriate tools
module load bioinfo-tools
module load MAFFT/7.407
module load trimAl/1.4.1
module load iqtree/2.0-rc2-omp-mpi
module load python/3.7.2
module load biopython/1.76-py3

#Command for running first part of bash pipeline (p1)
sh run_bash_p1.sh 
