#!/bin/bash

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/202305_micetimeline_48samples-MW-HdR/run_parameters_mice48_genome-n-ERCC.sh




################################################################################

# The following will -- given the settings in the parameter config files given above 
# -- map the datasets simulatenously to the ERCC spike ins and genome.

cd /hpc/hub_oudenaarden/mwehrens/fastq/202305_micetimeline_48samples-MW-HdR

mv Sample-Hesther-de-Ruiter_AACLJ2GM5_S18_L001_I1_001.fastq.gz Sample-Hesther-de-Ruiter_AACLJ2GM5_S18_L001_I1.fastq.gz
mv Sample-Hesther-de-Ruiter_AACLJ2GM5_S18_L001_R1_001.fastq.gz Sample-Hesther-de-Ruiter_AACLJ2GM5_S18_L001_R1.fastq.gz
mv Sample-Hesther-de-Ruiter_AACLJ2GM5_S18_L001_R2_001.fastq.gz Sample-Hesther-de-Ruiter_AACLJ2GM5_S18_L001_R2.fastq.gz

mkdir mapping

# S1
lib=Sample-Hesther-de-Ruiter_AACLJ2GM5_S18_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=10G
jobname=mice48
sbatch --job-name=${jobname} --time=72:00:00   --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh







#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 