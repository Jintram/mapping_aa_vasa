#!/bin/bash

# For dataset WANG2 N2, script failed because it ran out of 
# hard disk space ..
# So continue the job with 

# See
# https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/FirstTimeUsers#Submission_of_jobs
# for an overview of specs

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/run_parameters_Wang_N2-pt2.sh

lib=GSM2970358_N2_LV_cat
sbatch --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:1999999M --time=14-00:00:00  --mem=383G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh
  # note that there's only 3 nodes with this much space on the temp drive ..




#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 