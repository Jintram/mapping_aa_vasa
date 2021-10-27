#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/run_parameters_Wang.sh

lib=GSM3449620_N14_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=GSM2970358_N2_LV_cat
sbatch --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=GSM2970362_N3_LV_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=GSM2970366_N4_LV_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=GSM3449619_N13_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=GSM2970359_N6_LA_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=GSM2970363_N10_LA_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh


echo "omitted 2nd half"
if [[ "true" == "false" ]]; then
  lib=GSM2970367_N12_LA_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

  lib=GSM3449620_N14_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

  lib=GSM2970360_N5_LV_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

  lib=GSM2970364_N11_LA_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

  lib=GSM2970368_N8_LA_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

  lib=GSM2970361_N1_LV_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

  lib=GSM2970365_N9_LA_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

  lib=GSM2970369_N7_LA_cat
  sbatch  --output=slurm-${lib}-%x.%j.out --job-name=${lib} -c 8 --gres=tmpspace:999999M --time=14-00:00:00  --mem=255G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh
fi






#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 