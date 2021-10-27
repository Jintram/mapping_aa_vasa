#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give input (1) root file; (2) protocol [celseq1, celseq2, nla, scarsc]; (3) reference [mouse]"
    exit
fi

p2s=/hpc/hub_oudenaarden/aalemany/bin/mapandgo2
source /hpc/hub_oudenaarden/aalemany/virtualEnvironments/venv36/bin/activate

in=$1
out=${in%_*_S*_L*}
protocol=$2
reference=$3

# merge data
#${p2s}/mergeLanes.sh $in $out

# extract barcodes
#${p2s}/extractBC.sh $out ${protocol}

# trim
#${p2s}/trim.sh ${out}_cbc.fastq.gz

# map with bwa
#${p2s}/mapbwa.sh ${out}_cbc_trimmed.fq.gz ${out}_cbc_trimmed_bwa $reference

# map with star
${p2s}/mapstar.sh ${out}_cbc_trimmed.fq.gz ${out}_cbc_trimmed_star $reference

# create cout tables from bwa map 
#python ${p2s}/tablator_bwa.py ${out}_cbc_trimmed_bwa.bam

# create count tables from star map
#python ${p2s}/tablator_star.py ${out}_cbc_trimmed_starAligned.toTranscriptome.out.bam
#${p2s}/RNAvel_tables.sh ${out}_cbc_trimmed_starAligned.sortedByCoord.out.bam ${reference} ${out}_cbc_trimmed_star



