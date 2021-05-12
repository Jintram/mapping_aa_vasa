#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/run_parameters_all.sh

lib=HUB-JE-010_HGVN3BGX9_S1_cat
echo "skipping HUB-JE-010_HGVN3BGX9_S1_cat since already done (re-activate later!)"
#sbatch --job-name=derHammer --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=HUB-AL-s001_HG25TBGXF_S5_cat
sbatch --job-name=AL1 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=HUB-AL-s002_HG25TBGXF_S6_cat
sbatch --job-name=AL2 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE1_AHVJLVBGX3_S9_cat_
sbatch --job-name=JE1 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=HUB-JE-011_HGVN3BGX9_S2_cat
sbatch --job-name=JE11 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=rJE1_AHFL7NBGX5_S3_cat
sbatch --job-name=rJE1 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE2_AHY3WGBGX3_S1_cat
sbatch --job-name=JE2 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE3_AHY3WGBGX3_S2_cat
sbatch --job-name=JE3 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE4_AHFL7NBGX5_S4_cat
sbatch --job-name=JE4 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE5_AHFL77BGX5_S6_cat
sbatch --job-name=JE5 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE6_AHFL77BGX5_S7_cat
sbatch --job-name=JE6 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE7_AHFL7NBGX5_cat
sbatch --job-name=JE7 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=JE8_AHFL7NBGX5_S17_cat
sbatch --job-name=JE8 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=HUB-MW-005_AH32W2BGX9_S5_cat
sbatch --job-name=MW5 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=HUB-MW-006_AH32W2BGX9_S6_cat
sbatch --job-name=MW6 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=HUB-MW-007_HC3GFBGX9_S6_cat
sbatch --job-name=MW7 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

lib=HUB-MW-008_HC3GFBGX9_S7_cat
sbatch --job-name=MW8 --time=48:00:00  --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh


#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 