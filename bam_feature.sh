#!/bin/bash

#SBATCH -J gtex_download
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nickpowe@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=128G

export PATH=/filepath/subread/subread-2.0.1-source/bin/:$PATH

featureCounts -J -T 4 -s 0 -a /filepath/GTEx/gtf/gencode.v37.annotation.gtf.gz -o /filepath/GTEx/feature_counts_rna_seq_liver/liver_read_matrix2.txt /filepath/GTEx/bam_rna_seq_liver/*bam

featureCounts -T 4 -s 0 -a /filepath/GTEx/gtf/gencode.v37.annotation.gtf.gz -o /filepath/GTEx/feature_counts_rna_seq_sigmoid/sigmoid_read_matrix.txt /filepath/GTEx/bam_rna_seq_sigmoid/*bam

featureCounts -T 4 -s 0 -a /filepath/GTEx/gtf/gencode.v37.annotation.gtf.gz -o /filepath/GTEx/feature_counts_rna_seq_cerebellum/cerebellum_read_matrix.txt /filepath/GTEx/bam_rna_seq_cerebellum/*bam