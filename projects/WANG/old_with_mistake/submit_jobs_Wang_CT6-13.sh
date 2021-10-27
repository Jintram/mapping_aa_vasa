#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/run_parameters_Wang_CT6-13.sh

lib=GSM3449619_N13_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} --time=96:00:00  --mem=100G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=GSM2970359_N6_LA_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} --time=96:00:00  --mem=100G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh


#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 