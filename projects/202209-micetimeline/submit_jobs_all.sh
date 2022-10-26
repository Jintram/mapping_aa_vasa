#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/202209_micetimeline/run_parameters_micetest.sh



lib=Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=1G
sbatch --job-name=tinymice --time=48:00:00   --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=10G
sbatch --job-name=tinymice --time=48:00:00   --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh


################################################################################

# I later realized we also want to include the ERCC spike in reads.
# At a later point in time, I need to just merge this with the chromosome..

################################################################################

# Note that this file overwrites the genome mapping files!
# --> I need to update the ref genome such that they contain both 
#     mice genome and the ERCC reads.. (or some other solution)

################################################################################

# Repeat for ERCC reads

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/202209_micetimeline/run_parameters_ERCC.sh

lib=Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=1G
sbatch --job-name=tinymice --time=48:00:00  --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh
# dependency=15000394
# sbatch --dependency=afterok:${dependency} --job-name=tinymice --time=48:00:00  --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh 

lib=Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=10G
sbatch --job-name=tinymice --time=48:00:00  --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh
# dependency=15000394
# sbatch  --dependency=afterok:${dependency} --job-name=tinymice --time=48:00:00  --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh 

################################################################################




#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 