#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/WANG4/run_parameters_Wang.sh

# can already run this set, since it's quite small
lib=p.N2.plate.97493.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N1.plate.97474.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N1.plate.97474.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N2.plate.97452.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N2.plate.97452.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N2.plate.97493.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N2.plate.97493.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N3.plate.97438.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N3.plate.97438.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N4.plate.97461.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N4.plate.97461.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N5.plate.97458.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N5.plate.97458.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N13.plate.100355.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N13.plate.100355.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N14.plate.104720.part.1_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=p.N14.plate.104720.part.2_cat
sbatch  --output=slurm-${lib}-%x.%j.out --job-name=map_${lib} -c 8 --gres=tmpspace:499999M --time=14-00:00:00  --mem=128G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh



# currently still running:
# p.N5.plate.97458.part.2  8213570



#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 