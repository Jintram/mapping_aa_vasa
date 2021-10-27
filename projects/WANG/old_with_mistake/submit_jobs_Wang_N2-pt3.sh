#!/bin/bash

# For dataset WANG2 N2, script failed because it ran out of 
# hard disk space ..
# Also when I increased resources. So I now split the job up according to 
# barcodes, in 10 ±equal parts.

# See
# https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/FirstTimeUsers#Submission_of_jobs
# for an overview of specs

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/run_parameters_Wang_N2-pt3.sh

for thepart in part.0. part.1. part.2. part.3. part.4. part.5. part.6. part.7. part.8. part.9.
do
  lib=${thepart}GSM2970358_N2_LV_cat
  sbatch --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 1 --gres=tmpspace:999999M --time=14-00:00:00  --mem=150G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh
done





#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 