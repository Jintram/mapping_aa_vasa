#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
outputdir=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/homeMadeBed/
outputfile=homeMade_IntronsExons__Homo_sapiens.GRCh38.93_
output=${outputdir}${outputfile}

mkdir $outputdir

sbatch --job-name=BedsBurning --time=48:00:00  --mem=30G --export=ALL,gtffile="${gtffile}",output="${output}" ${script_path}/setting_up/convert_gtf_to_bed.sh

#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 
