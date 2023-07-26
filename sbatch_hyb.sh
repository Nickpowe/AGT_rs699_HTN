#!/bin/bash

#SBATCH -J hyb
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nickpowe@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=32G


module load bowtie2

cd /filepath/hyb

hyb analyse in=/filepath/name_of_file.fastq.gz db=hyb_ref_2020

