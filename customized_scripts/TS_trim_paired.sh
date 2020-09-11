#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give (1) input fastq file; (3) path2trimgalore; (4) path2cutadapt;"
  exit
fi

file2trim=$1
path2trimgalore=$2
path2cutadapt=$3
cutoff=10 # cutoff length to remove reads

# trim adaptors
${path2trimgalore}/trim_galore --paired --path_to_cutadapt ${path2cutadapt}/cutadapt ${file2trim}_R1_cbc.fastq.gz ${file2trim}_R2_cbc.fastq.gz --length $cutoff --retain_unpaired
mv ${file2trim}_R1_cbc.fastq.gz_trimming_report.txt ${file2trim}_R1_trimming_report.txt
mv ${file2trim}_R2_cbc.fastq.gz_trimming_report.txt ${file2trim}_R2_trimming_report.txt 
  # Note1: "For adapter trimming, Trim Galore! uses the first 13 bp of Illumina 
  # standard adapters ('AGATCGGAAGAGC') by default (suitable for both ends of 
  # paired-end libraries), but accepts other adapter sequence, too"
  # (trim galore documentation at https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  # See also ./Docs/Trim_Galore_User_Guide.md in trim galore directory;
  # Note2: default behavior is to remove pairs of which one read is <20nt,
  # this is compatible with targeted reads, where the the biological read 
  # will be ±41nt (75-20-14).
  # HOWEVER: with a cutoff of 20nt, ±10% of reads get removed, which seems
  # very high; I now lowered cutoff length to 20nt, and also used the retain_unpaired
  # option.
  # It is however unclear to me why some of these reads get removed when I look them up,
  # it seems both reads are >20nt, not of poor quality, and not containing ill primer..
  # --> TODO --> look into this..

# trim homopolymers 
${path2trimgalore}/trim_galore --paired --path_to_cutadapt ${path2cutadapt}/cutadapt -a AAAAAAAAAAAAAAAAAAAA --no_report_file ${file2trim}_R1_cbc_val_1.fq.gz ${file2trim}_R2_cbc_val_2.fq.gz --length $cutoff --retain_unpaired
mv ${file2trim}_R1_cbc_val_1_val_1.fq.gz ${file2trim}_R1_cbc_val_1_HA.fq.gz
mv ${file2trim}_R2_cbc_val_2_val_2.fq.gz ${file2trim}_R2_cbc_val_2_HA.fq.gz

${path2trimgalore}/trim_galore --paired --path_to_cutadapt ${path2cutadapt}/cutadapt -a TTTTTTTTTTTTTTTTTTTT --no_report_file ${file2trim}_R1_cbc_val_1_HA.fq.gz ${file2trim}_R2_cbc_val_2_HA.fq.gz --length $cutoff --retain_unpaired
mv ${file2trim}_R1_cbc_val_1_HA_val_1.fq.gz ${file2trim}_R1_cbc_val_1_HAT.fq.gz
rm ${file2trim}_R1_cbc_val_1_HA_val_1.fq.gz
mv ${file2trim}_R2_cbc_val_2_HA_val_2.fq.gz ${file2trim}_R2_cbc_val_2_HAT.fq.gz
rm ${file2trim}_R2_cbc_val_2_HA_val_2.fq.gz

${path2trimgalore}/trim_galore --paired --path_to_cutadapt ${path2cutadapt}/cutadapt -a CCCCCCCCCCCCCCCCCCCC --no_report_file ${file2trim}_R1_cbc_val_1_HAT.fq.gz ${file2trim}_R2_cbc_val_2_HAT.fq.gz --length $cutoff --retain_unpaired
mv ${file2trim}_R1_cbc_val_1_HAT_val_1.fq.gz ${file2trim}_R1_cbc_val_1_HATC.fq.gz
rm ${file2trim}_R1_cbc_val_1_HAT_val_1.fq.gz
mv ${file2trim}_R2_cbc_val_2_HAT_val_2.fq.gz ${file2trim}_R2_cbc_val_2_HATC.fq.gz
rm ${file2trim}_R2_cbc_val_2_HAT_val_2.fq.gz

${path2trimgalore}/trim_galore --paired --path_to_cutadapt ${path2cutadapt}/cutadapt -a GGGGGGGGGGGGGGGGGGGG --no_report_file ${file2trim}_R1_cbc_val_1_HATC.fq.gz ${file2trim}_R2_cbc_val_2_HATC.fq.gz --length $cutoff --retain_unpaired
mv ${file2trim}_R1_cbc_val_1_HATC_val_1.fq.gz ${file2trim}_R1_cbc_val_1_HATCG.fq.gz
rm ${file2trim}_R1_cbc_val_1_HATC_val_1.fq.gz
mv ${file2trim}_R2_cbc_val_2_HATC_val_2.fq.gz ${file2trim}_R2_cbc_val_2_HATCG.fq.gz
rm ${file2trim}_R2_cbc_val_2_HATC_val_2.fq.gz





#java -jar ${path2trimmomatic}/trimmomatic-0.36.jar SE -phred33 ${file2trim%.fastq.gz}_trimmed.fq.gz ${file2trim%.fastq.gz}_trimHomo.fq ILLUMINACLIP:/hpc/hub_oudenaarden/fsalmen/barcodes/homopolymers.fa:0::0:0 MINLEN:20

#gzip ${file2trim%.fastq.gz}_trimHomo.fq


