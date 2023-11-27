#!/bin/bash

# See also https://github.com/Jintram/mapping_aa_vasa for a guide to these scripts.


script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/202308_HFTFs_mapping/run_parameters_HFTFs_genome-n-ERCC.sh




################################################################################

# The following will -- given the settings in the parameter config files given above 
# -- map the datasets simulatenously to the ERCC spike ins and genome.

cd /hpc/hub_oudenaarden/mwehrens/fastq/202308_HFTFs_mapping/

mv i15-MW_AACTH7HM5_S12_L001_I1_001.fastq.gz i15-MW_AACTH7HM5_S12_L001_I1.fastq.gz
mv i15-MW_AACTH7HM5_S12_L001_R1_001.fastq.gz i15-MW_AACTH7HM5_S12_L001_R1.fastq.gz
mv i15-MW_AACTH7HM5_S12_L001_R2_001.fastq.gz i15-MW_AACTH7HM5_S12_L001_R2.fastq.gz

mv i16-MW_AACTH7HM5_S13_L001_I1_001.fastq.gz i16-MW_AACTH7HM5_S13_L001_I1.fastq.gz
mv i16-MW_AACTH7HM5_S13_L001_R1_001.fastq.gz i16-MW_AACTH7HM5_S13_L001_R1.fastq.gz
mv i16-MW_AACTH7HM5_S13_L001_R2_001.fastq.gz i16-MW_AACTH7HM5_S13_L001_R2.fastq.gz

mkdir mapping

# idx 15
lib=i15-MW_AACTH7HM5_S12_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=10G
jobname=HFTF15
sbatch --job-name=${jobname} --time=72:00:00   --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

# idx 16
lib=i16-MW_AACTH7HM5_S13_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=10G
jobname=HFTF16
sbatch --job-name=${jobname} --time=72:00:00   --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh


################################################################################

# Overexpression of our genes of interest was performed using a plasmid construct 
# that results in an mRNA molecule with an artificial 3' end. So to measure
# expression of the over expressed gene from that plasmid, let's also
# check mapping to that plasmid.

# First "install" the plasmid sequence; using custom .fa and .gtf files, that I made.
# local
cd /Users/m.wehrens/Data/2023_08_HFTFs/ref_sequences_custom
scp pLVX-IRES-Hyg.* mwehrens@gw2hpcs05:/hpc/hub_oudenaarden/mwehrens/ref/custom/pLVX-IRES-Hyg
# @HPC
# See file ./mapping_aa_private_vasa/STAR/submit_index_star_pLVX-IRES-Hyg.sh
# Using, ./mapping_aa_private_vasa/STAR/generate_index_STAR.sh

# Upload adjusted files
scp /Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/projects/202308-HFTFs_Alejandro-Martijn/run_parameters_HFTFs_pLVX-IRES-Hyg.sh mwehrens@gw2hpcs05:/hpc/hub_oudenaarden/mwehrens/fastq/202308_HFTFs_mapping/pLVX-IRES-Hyg

# Now map to this
# Shared with above
script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
# specific for pLVX-IRES-Hyg
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/202308_HFTFs_mapping/pLVX-IRES-Hyg/run_parameters_HFTFs_pLVX-IRES-Hyg.sh

# idx 15
lib=i15-MW_AACTH7HM5_S12_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=10G
jobname=HFTF15
sbatch --job-name=${jobname} --time=72:00:00   --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

# idx 16
lib=i16-MW_AACTH7HM5_S13_L001 # Fastq files name, without _Rx.fastq.gz
tmpspace=10G
jobname=HFTF16
sbatch --job-name=${jobname} --time=72:00:00   --gres=tmpspace:${tmpspace} --mem=50G --export=ALL,general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

# Then see
# /Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/projects/202308-HFTFs_Alejandro-Martijn/submit_Counttables_umiCounts_HFTFs_pLVX-IRES-Hyg_filteredanddiscardmulti.sh



#sbatch --time=05:15:00  --mem=10G L_TS_submit_vasaplate_map_LOCAL.sh 


# Now copy them back
cd 
scp mwehrens@gw2hpcs05:"/hpc/hub_oudenaarden/mwehrens/fastq/202308_HFTFs_mapping/pLVX-IRES-Hyg/*counts.tsv.gz" /Users/m.wehrens/Data/2023_08_HFTFs/count-tables/mapped-pLVX-IRES-Hyg
scp mwehrens@gw2hpcs05:"/hpc/hub_oudenaarden/mwehrens/fastq/202308_HFTFs_mapping/pLVX-IRES-Hyg/*assigned_sorted_MW.bam*" /Users/m.wehrens/Data/2023_08_HFTFs/count-tables/mapped-pLVX-IRES-Hyg




















