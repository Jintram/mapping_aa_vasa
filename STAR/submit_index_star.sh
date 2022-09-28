#!/bin/bash

fafile=./ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtffile=./ensembl/Homo_sapiens.GRCh38.93.gtf
overhang=74

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 

fafile=./ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtffile=./ensembl/Homo_sapiens.GRCh38.93.gtf
overhang=49

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 



