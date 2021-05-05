#!/bin/bash

currentpath=/Volumes/workdrive_m.wehrens_hubrecht/Fastq__raw_files/HCM_SCS/2021_combined_cat_txt/

for current_file in JE2_AHY3WGBGX3_S1 JE3_AHY3WGBGX3_S2 JE4_AHFL7NBGX5_S4 JE7_AHFL7NBGX5 JE8_AHFL7NBGX5_S17 HUB-JE-010_HGVN3BGX9_S1 HUB-JE-011_HGVN3BGX9_S2 HUB-MW-005_AH32W2BGX9_S5 HUB-MW-006_AH32W2BGX9_S6 HUB-MW-007_HC3GFBGX9_S6 HUB-MW-008_HC3GFBGX9_S7 JE1_AHVJLVBGX3_S9 JE5_AHFL77BGX5_S6 JE6_AHFL77BGX5_S7 HUB-AL-s001_HG25TBGXF_S5 HUB-AL-s002_HG25TBGXF_S6 rJE1_AHFL7NBGX5_S3
do
  echo processing ${current_file}
	cat ${currentpath}*${current_file}*R1*.fastq | wc -l
  cat ${currentpath}*${current_file}*R2*.fastq | wc -l
done