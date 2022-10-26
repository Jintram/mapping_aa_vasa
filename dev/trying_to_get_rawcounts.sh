#!/bin/sh

f=Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT.nonRibo_E99_Aligned.out.bam.featureCounts
samtools view -h -o $f.sam $f.bam
f12=Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT.nonRibo_E99_Aligned.out.bam.featureCounts
samtools view -h -o ${f12}.sam ${f12}.bam

head -1000 $f.sam > head_$f.sam

# (Random awk commands; forget about this)
#cat head_$f.sam | awk '{print $1}'
#cat head_$f.sam | awk '{sub(/.*-/,"",$2);print $2,$4}' 

# The following command doesn't demultiplex BCs
featureCounts -a ${gtffile} -o raw-counts -s 1 $bamfile.bam

# byReadGroup seems to be convenient
f2=Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT.nonRibo_E99_Aligned.out
samtools view -h -o $f2.sam $f2.bam




# The "ignore-umi" options doesn't lead to reaw counts; just collapses everything --> all reads are 1
umi_tools count --ignore-umi --cell-tag=SM --extract-umi-method=tag --umi-tag=RX --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${bamfile}.assigned_sorted_MW.bam -S ${bamfile}.counts-raw.tsv.gz

