#!/bin/bash

#SBATCH -t 02:00 		# hours:minutes runlimit after which job will be killed
#SBATCH -c 6 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 16G
#SBATCH --job-name STAR_index 		# Job name

# Following tutorial from
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

#cd /n/scratch2/username/

#module load gcc/6.2.0 star/2.5.2b

gunzip *

/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin/STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ./STAR-indexed-74/ \
--genomeFastaFiles ./Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--sjdbGTFfile ./Homo_sapiens.GRCh38.81.gtf.gz \
--sjdbOverhang 74
