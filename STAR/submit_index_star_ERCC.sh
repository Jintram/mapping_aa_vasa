#!/bin/bash

# Mouse genome for test run of RNA library seq 

fafile=./ERCC92.fa
gtffile=./ERCC92.gtf
overhang=150

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 