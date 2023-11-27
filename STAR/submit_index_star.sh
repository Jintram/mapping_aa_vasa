#!/bin/bash

###

fafile=./ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtffile=./ensembl/Homo_sapiens.GRCh38.93.gtf
overhang=74

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 

###

fafile=./ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtffile=./ensembl/Homo_sapiens.GRCh38.93.gtf
overhang=49

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 

###

# cd /hpc/hub_oudenaarden/mwehrens/ref/GRCh38.110
fafile=./Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gtffile=./Homo_sapiens.GRCh38.110.gtf.gz
overhang=100

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 

###

# Add ERCC to genome
cd /hpc/hub_oudenaarden/mwehrens/ref/
zcat /hpc/hub_oudenaarden/mwehrens/ref/ERCC/ERCC92.fa.gz ./GRCh38.110/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ./GRCh38.110.ERCC/Homo_sapiens.GRCh38.110.dna.primary_assembly_plus_ERCC.fa
zcat /hpc/hub_oudenaarden/mwehrens/ref/ERCC/ERCC92.gtf.gz ./GRCh38.110/Homo_sapiens.GRCh38.110.gtf.gz > ./GRCh38.110.ERCC/Homo_sapiens.GRCh38.110_plus_ERCC.gtf
gzip ./GRCh38.110.ERCC/Homo_sapiens.GRCh38.110.dna.primary_assembly_plus_ERCC.fa
gzip ./GRCh38.110.ERCC/Homo_sapiens.GRCh38.110_plus_ERCC.gtf

cd GRCh38.110.ERCC
fafile=./Homo_sapiens.GRCh38.110.dna.primary_assembly_plus_ERCC.fa
gtffile=./Homo_sapiens.GRCh38.110_plus_ERCC.gtf
overhang=100

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 
