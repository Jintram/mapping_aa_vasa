#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/run_parameters.sh

lib=HUB-JE-010_HGVN3BGX9_S1_cat
sbatch --job-name=derHammer --time=24:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 