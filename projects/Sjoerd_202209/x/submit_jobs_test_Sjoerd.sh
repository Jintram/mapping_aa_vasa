#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/Sjoerd/run_parameters_Sjoerd.sh


lib=head_clone-12-no-dox-rep1_AAC3KKWM5_S19_L001
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=test --time=05:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 