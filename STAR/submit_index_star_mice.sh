#!/bin/bash

# Mouse genome for test run of RNA library seq 

fafile=./ensembl/Mus_musculus.GRCm39.dna.primary_assembly.fa
gtffile=./ensembl/Mus_musculus.GRCm39.107.gtf
overhang=150

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 