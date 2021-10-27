#!/bin/bash

my_directory=/Volumes/fastq_m.wehrens/Mapping/WANG4/
file=p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out

# samtools sort -T /tmp/aln.sorted -o aln.sorted.bam aln.bam
samtools sort -o ${my_directory}${file}.sorted.bam ${my_directory}${file}.bam
samtools index ${my_directory}${file}.sorted.bam