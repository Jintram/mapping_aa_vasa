#!/bin/bash
# Part 2 of the mapping pipeline: 
# BWA mapping to Ribosomal reference reads. 
# Produces  an output fastq.gz file containing non-ribosomal reads;
# and a bam file containing ribosomal reads. This last file contains also
# reads that map with very low qualities. 

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "Please, give (1) general param file, (2) run param file, (3) fq file to map, (4) out-prefix"    
    exit
fi

general_parameter_filepath=$1
run_parameter_filepath=$2
fq=$3
out=$4

source $general_parameter_filepath
source $run_parameter_filepath
current_dir=$(pwd)
echo "Going to $outdir"
cd $outdir

# mapping short reads
echo "mapping rRNA, short"
${p2bwa}/bwa aln ${riboref} ${fq} > aln_${fq%.f*q}.sai 
${p2bwa}/bwa samse  ${riboref}  aln_${fq%.f*q}.sai ${fq} > aln-temp.out 
${p2samtools}/samtools view -Sb aln-temp.out  > ${out}.aln-ribo.bam 
rm aln-temp.out # i don't know why, but locally it only worked with creating an intermediate file

# mapping normal reads
echo "mapping rRNA, long"
${p2bwa}/bwa mem -t 8 -h 15 ${riboref} ${fq} > mem-temp.out
${p2samtools}/samtools view -Sb mem-temp.out > ${out}.mem-ribo.bam 
rm mem-temp.out

echo "merging two maps" 
${p2samtools}/samtools merge -f -n -r -h ${out}.aln-ribo.bam ${out}.all-ribo.bam ${out}.aln-ribo.bam ${out}.mem-ribo.bam 
${p2samtools}/samtools view -H ${out}.aln-ribo.bam

rm ${out}.aln-ribo.bam ${out}.mem-ribo.bam aln_${fq%.f*q}.sai

echo "sorting"
${p2samtools}/samtools sort -n -O bam -T temp_${out} -o ${out}.nsorted.all-ribo.bam ${out}.all-ribo.bam
rm ${out}.all-ribo.bam

echo "riboread selection"
#pythonbin=/Users/m.wehrens/anaconda3/bin/python
$pythonbin ${p2s}/TS_riboread-selection.py ${out}.nsorted.all-ribo.bam $stranded ${out}

echo "end of script"
exit

#if [ $stranded == "n" ]
#then
#    # select reads that don't map and create a new fastq file
#    ${p2samtools}/samtools view -f 4 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq & 
#    # select reads that map (even low qualities, and create a new bam file
#    ${p2samtools}/samtools view -F 4 -Sb ${out}.ribo.sam > ${out}.ribo.bam & 
#elif [ $stranded == "y" ]
#then
#    ${p2samtools}/samtools view -f 16 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq
#    ${p2samtools}/samtools view -f 4 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq &
#    ${p2samtools}/samtools view -f 0 -F 4 -Sb ${out}.ribo.sam > ${out}.ribo.bam &
#fi

#wait

# zip fastq file and delete sam file from first mapping
#gzip ${out}.nonRibo.fastq &
#rm ${out}.ribo.sam

#wait

