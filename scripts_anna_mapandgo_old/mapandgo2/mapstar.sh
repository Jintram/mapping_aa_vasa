#!/bin/bash

#path2star=/hpc/hub_oudenaarden/avo/nascent/STAR-2.5.3a/bin/Linux_x86_64

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) fastq file to map"
    echo "2) root for output file (no .sam or .bam extension)"
    echo "3) reference genome (mouse, human, zebrafish, spiny, celegans)"
    echo "4) path to STAR"
    exit
fi

file2map=$1
outfq=$2
ref=$3
path2star=$4

if [ $ref == 'mouse' ]
then
    ref=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/star_index_60
elif [ $ref == 'human' ]
then
    ref=/hpc/hub_oudenaarden/group_references/ensembl/93/homo_sapiens/star_index_60
elif [ $ref == 'zebrafish' ]
then
    ref=/hpc/hub_oudenaarden/group_references/ensembl/93/danio_rerio/star_index_75
elif [ $ref == 'celegans' ]
then
    ref=hpc/hub_oudenaarden/group_references/ensembl/94/c_elegans/star_index_75
elif [ $ref == 'spiny' ]
then
    ref=/hpc/hub_oudenaarden/aalemany/henriette/ref_genome_spiny/star_index_75
fi

#${path2star}/STAR --runThreadN 12 --genomeDir $ref --readFilesIn ${file2map} --readFilesCommand zcat --outFileNamePrefix ${outfq} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM
#rm -r ${outfq}_STARtmp

${path2star}/STAR --runThreadN 12 --genomeDir $ref --readFilesIn ${file2map} --readFilesCommand zcat --outFileNamePrefix ${outfq} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif --outFilterMultimapNmax 1 
rm -r ${outfq}_STARtmp


