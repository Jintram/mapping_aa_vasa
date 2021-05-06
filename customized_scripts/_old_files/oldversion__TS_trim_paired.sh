#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give (1) input fastq file; (2) path2trimgalore; (3) path2cutadapt;"
  exit
fi

file2trim=$1
path2trimgalore=$2
path2cutadapt=$3

# trim adaptors
${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt ${file2trim}
mv ${file2trim}_trimming_report.txt ${file2trim%.fastq.gz}_trimming_report.txt
  # note that "For adapter trimming, Trim Galore! uses the first 13 bp of Illumina 
  # standard adapters ('AGATCGGAAGAGC') by default (suitable for both ends of 
  # paired-end libraries), but accepts other adapter sequence, too"
  # (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  # Note the default behavior is to discard paired reads when one of the reads
  # is <20nt. This is suitable for the targeted sequencing also, as the biological 
  # region is 42nt.

# trim homopolymers 
${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a AAAAAAAAAAAAAAAAAAAA --no_report_file ${file2trim%.fastq.gz}_trimmed.fq.gz
mv ${file2trim%.fastq.gz}_trimmed_trimmed.fq.gz ${file2trim%.fastq.gz}_trimmed_HA.fq.gz

${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a TTTTTTTTTTTTTTTTTTTT --no_report_file ${file2trim%.fastq.gz}_trimmed_HA.fq.gz
mv ${file2trim%.fastq.gz}_trimmed_HA_trimmed.fq.gz ${file2trim%.fastq.gz}_trimmed_HAT.fq.gz
rm ${file2trim%.fastq.gz}_trimmed_HA.fq.gz

${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a CCCCCCCCCCCCCCCCCCCC --no_report_file ${file2trim%.fastq.gz}_trimmed_HAT.fq.gz
mv ${file2trim%.fastq.gz}_trimmed_HAT_trimmed.fq.gz ${file2trim%.fastq.gz}_trimmed_HATC.fq.gz
rm ${file2trim%.fastq.gz}_trimmed_HAT.fq.gz

${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a GGGGGGGGGGGGGGGGGGGG --no_report_file ${file2trim%.fastq.gz}_trimmed_HATC.fq.gz
mv ${file2trim%.fastq.gz}_trimmed_HATC_trimmed.fq.gz ${file2trim%.fastq.gz}_trimmed_homoATCG.fq.gz
rm ${file2trim%.fastq.gz}_trimmed_HATC.fq.gz





#java -jar ${path2trimmomatic}/trimmomatic-0.36.jar SE -phred33 ${file2trim%.fastq.gz}_trimmed.fq.gz ${file2trim%.fastq.gz}_trimHomo.fq ILLUMINACLIP:/hpc/hub_oudenaarden/fsalmen/barcodes/homopolymers.fa:0::0:0 MINLEN:20

#gzip ${file2trim%.fastq.gz}_trimHomo.fq


