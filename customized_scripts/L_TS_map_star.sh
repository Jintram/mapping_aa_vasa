#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give (1) general param file, (2) run param file, (3) input fastq file (only one), (4) prefix for output file names"
    #echo "Please, give as input:"
    #echo "1) path to STAR"
    #echo "2) path to samtools"
    #echo "3) genome folder"
    #echo "4) input fastq file (only one)"
    #echo "5) prefix for output file names"
    exit
fi



################################################################################
general_parameter_filepath=$1
run_parameter_filepath=$2
inputfq=$3
outprefix=$4

source $general_parameter_filepath
source $run_parameter_filepath

current_dir=$(pwd)

if [[ $TMPDIR == "" ]]; then
  echo "Going to $outdir"
  cd $outdir
else
  echo "Going to $TMPDIR"
  cd $TMPDIR
fi

################################################################################

# For local use:
# - zcat doesn't work for some reason on local mac
# - not multi-threaded
#gunzip ${inputfq}.gz
#${p2star}/STAR --genomeDir ${genome} --readFilesIn ${inputfq} --readFilesCommand cat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix ${outprefix}

${p2star}/STAR --runThreadN 8 --genomeDir ${genome} --readFilesIn ${inputfq}.gz --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix ${outprefix}

rm -r ${outprefix}_STARtmp
rm ${outprefix}Log.progress.out

mv ${outprefix}Log.out ${outprefix}Log.txt
mv ${outprefix}Log.final.out ${outprefix}Log.final.txt

#${p2samtools}/samtools view -q 255 ${outprefix}Aligned.sortedByCoord.out.bam -b -o ${outprefix}Aligned.sortedByCoord.filtered.bam
#rm ${outprefix}Aligned.sortedByCoord.out.bam

################################################################################

if [[ -f ${outprefix}Aligned.out.bam ]]; then
  
  if [[ $nocleanup = "" ]]; then
    echo "cleaning up some files" # prevent this by setting "nocleanup"
    #rm ${inputfq}
    rm ${inputfq}.gz
  fi
  
else
  
  echo "output files not found, returning error"
  exit 1
  
fi

# copy output files from temp dir to output dir 
if [[ $TMPDIR != "" ]]; then
  
  cp ${outprefix}Aligned.out.bam $outdir/
  
fi

################################################################################






