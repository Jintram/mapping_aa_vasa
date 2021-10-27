#!/bin/bash


for i in {0..9}
do
    file=part.${i}.GSM2970358_N2_LV_cat_nc.nonRibo_E99_Aligned.out
    samtools view -S -b ${file}.sam > ${file}.bam &
done

echo "all jobs running, now waiting .."

wait

echo "Conversion complete."

# The opposite,
# to convert from bam to sam:
# samtools view -h -o GSM2970358_N2_LV_cat_nc.nonRibo_E99_Aligned.out.sam GSM2970358_N2_LV_cat_nc.nonRibo_E99_Aligned.out.bam
# samtools view -h -o out.sam in.bam
