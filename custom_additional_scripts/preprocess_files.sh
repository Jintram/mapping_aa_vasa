#!/bin/bash

currentpath=/Volumes/workdrive_m.wehrens_hubrecht/sequencing_raw_files/Allelic_Pilot3_cat/

gunzip ${currentpath}*.gz

for current_file in MW-TS-S1A-Hybr-Bst_HTMH2BGXF_S25 MW-TS-S1-Hybr-NB_HTMH2BGXF_S23 MW-TS-S2A-TargOnly-Bst_HTMH2BGXF_S26 MW-TS-S2-TargOnly-NB_HTMH2BGXF_S24 MW-TS-S4A-Vasa-Bst_HTMH2BGXF_S27 MW-TS-S4-Vasa-NB_HTMH2BGXF_S28
do
	cat ${currentpath}${current_file}*R1*.fastq > ${currentpath}${current_file}_R1.fastq
  gzip ${currentpath}${current_file}_R1.fastq
  cat ${currentpath}${current_file}*R2*.fastq > ${currentpath}${current_file}_R2.fastq
  gzip ${currentpath}${current_file}_R2.fastq
done

