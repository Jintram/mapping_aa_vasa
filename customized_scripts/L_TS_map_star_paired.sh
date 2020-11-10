#!/bin/bash

if [ $# -ne 6 ]
then
    echo "Please, give as input:"
    echo "1) path to STAR"
    echo "2) path to samtools"
    echo "3) genome folder"
    echo "4) input fastq file R1 (only one)"
    echo "5) input fastq file R2"
    echo "6) prefix for output file names"
    exit
fi

p2star=$1
p2samtools=$2
genome=$3
inputfq1=$4
inputfq2=$5
outprefix=$6

gunzip ${inputfq1}.gz ${inputfq2}.gz
${p2star}/star --genomeDir ${genome} --readFilesIn ${inputfq1} ${inputfq2} --readFilesCommand cat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix ${outprefix}

rm -r ${outprefix}_STARtmp
rm ${outprefix}Log.progress.out

mv ${outprefix}Log.out ${outprefix}Log.txt
mv ${outprefix}Log.final.out ${outprefix}Log.final.txt

#${p2samtools}/samtools view -q 255 ${outprefix}Aligned.sortedByCoord.out.bam -b -o ${outprefix}Aligned.sortedByCoord.filtered.bam
#rm ${outprefix}Aligned.sortedByCoord.out.bam


