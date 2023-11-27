#!/bin/bash

# Set files to index
fafile=./pLVX-IRES-Hyg.fa
gtffile=./pLVX-IRES-Hyg.gtf
overhang=100

# Submit
sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 