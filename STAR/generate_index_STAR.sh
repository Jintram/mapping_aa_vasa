#!/bin/bash

#SBATCH --time 24:00:00                 # hours:minutes:secs runlimit after which job will be killed
#SBATCH -c 6            # number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#earlier using  --mem 100G
#SBATCH --mem-per-cpu 30G
#SBATCH --job-name brightSTAR           # Job name

# Following tutorial from
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

#cd /n/scratch2/username/

#module load gcc/6.2.0 star/2.5.2b

#fafile=./Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#gtffile=./Homo_sapiens.GRCh38.81.gtf.gz

if [[ "${fafile}" == "" ]]; then
  fafile=./ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa
  gtffile=./ensembl/Homo_sapiens.GRCh38.93.gtf
  overhang=74
fi

gunzip ${fafile}.gz
gunzip ${gtffile}.gz

/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin/STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ./STAR-indexed-${overhang}/ \
--genomeFastaFiles ${fafile} \
--sjdbGTFfile ${gtffile} \
--sjdbOverhang ${overhang}

# gzip ./ensembl/*