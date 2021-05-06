

# THIS SCRIPT WAS EDITED TO PERFORM MAPPING OF PAIRS.
# 
# HOWEVER THIS IS STRICTLY NOT NECESSARY, AND IT IS A BIT INCONVENIENT TO DEAL
# WITH PAIRED READS IN THE FRAMEWORK OF THESE SCRIPTS, SPECIFICALLY, riboread-selection.py
# USES THE MERGED .BAM OUTPUT FILES, WHICH THEN CONTAIN BOTH MAPPING FROM 


#!/bin/bash
# Part 2 of the mapping pipeline: 
# BWA mapping to Ribosomal reference reads. 
# Produces  an output fastq.gz file containing non-ribosomal reads;
# and a bam file containing ribosomal reads. This last file contains also
# reads that map with very low qualities. 

if [ $# -ne 8 ]
then
    echo "Please, give:"
    echo "(1) reference fasta file (mouse, human or full path)"
    echo "(2) fastq file R1 to map"
    echo "(3) fastq file R2 to map"
    echo "(4) prefix output"
    echo "(5) Path to bwa software"
    echo "(6) Path to samtools"
    echo "(7) Stranded protocol? (y/n)"
    echo "(8) Path to riboread-selection.py script"
    exit
fi

ref=$1
fq_R1=$2
fq_R2=$3
out=$4
p2bwa=$5
p2samtools=$6
stranded=$7
p2s=$8

# mapping short reads
# (for paired reads)
${p2bwa}/bwa index ${ref} 
${p2bwa}/bwa aln ${ref} ${fq_R1} > aln_${fq_R1%.f*q.gz}.sai
${p2bwa}/bwa aln ${ref} ${fq_R2} > aln_${fq_R2%.f*q.gz}.sai

${p2bwa}/bwa sampe  ${ref}  aln_${fq_R1%.f*q.gz}.sai aln_${fq_R2%.f*q.gz}.sai ${fq_R1} ${fq_R2} | ${p2samtools}/samtools view -Sb > ${out}.aln-ribo.bam &
  # test locally:
  #${p2bwa}/bwa sampe  ${ref}  aln_${fq_R1%.f*q.gz}.sai aln_${fq_R2%.f*q.gz}.sai ${fq_R1} ${fq_R2} > bwa.temp 
  #${p2samtools}/samtools view -Sb bwa.temp > ${out}.aln-ribo.bam &
  #rm bwa.temp # i don't know why, but locally it only worked with creating an intermediate file

# mapping normal reads
# (for paired reads)
${p2bwa}/bwa mem -t 8 -h 15 ${ref} ${fq_R1} ${fq_R2} | ${p2samtools}/samtools view -Sb > ${out}.mem-ribo.bam & 
  # test locally:
  #${p2bwa}/bwa mem -t 8 -h 15 ${ref} ${fq_R1} ${fq_R2} > mem.temp
  #${p2samtools}/samtools view -Sb mem.temp > ${out}.mem-ribo.bam # & 
  #rm mem.temp

wait

# merge both short and normal read mapping into one file
# Server command:
${p2samtools}/samtools merge -n -r -h ${out}.aln-ribo.bam --threads 8 ${out}.all-ribo.bam ${out}.aln-ribo.bam ${out}.mem-ribo.bam 
# For local testing:
# ${p2samtools}/samtools merge -f -n -r -h ${out}.aln-ribo.bam ${out}.all-ribo.bam ${out}.aln-ribo.bam ${out}.mem-ribo.bam 
  # use: 
  # p2samtools=/Users/m.wehrens/Software_custom/samtools-1.2
  # and have removed the "--threads 8" option
  # -f was added because multiple tests might require overwriting earlier files
rm ${out}.aln-ribo.bam ${out}.mem-ribo.bam aln_${fq%.f*q}.sai
${p2samtools}/samtools sort -n --threads 8 ${out}.all-ribo.bam -O BAM -o ${out}.nsorted.all-ribo.bam
# For local testing:
# ${p2samtools}/samtools sort -n -O bam -T temp_${out} -o ${out}.nsorted.all-ribo.bam  ${out}.all-ribo.bam
  # with
  # p2samtools=/Users/m.wehrens/Software_custom/samtools-1.2
rm ${out}.all-ribo.bam

${p2s}/TS_riboread-selection_paired.py ${out}.nsorted.all-ribo.bam $stranded ${out}
  # for local testing:
  # pythonbin=/Users/m.wehrens/anaconda3/bin/python
  # $pythonbin ${p2s}/TS_riboread-selection_paired.py ${out}.nsorted.all-ribo.bam $stranded ${out}
  #
  # For debugging purposes, it can be convenient to index the bam file,
  # to do that however, it needs to be sorted by position first
  # ${p2samtools}/samtools sort -O bam -T temp_${out} -o ${out}.isorted.all-ribo.bam ${out}.all-ribo.bam
  # ${p2samtools}/samtools index ${out}.isorted.all-ribo.bam
  # 
  # also convenient to convert the bam file to a sam file just to inspect the contents
  # samtools view -h -o MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001.Ribo.sam MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001.Ribo.bam
  #
  # Or convert sorted file to fastq file:
  # /Users/m.wehrens/Software_custom/bedtools2/bin/bamToFastq -i MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001.isorted.all-ribo.bam -fq MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001.isorted.all-ribo.fq


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

