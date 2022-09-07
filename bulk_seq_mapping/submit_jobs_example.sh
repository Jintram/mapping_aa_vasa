#!/bin/bash

################################################################################

# It's most convenient to simply paste the contents of this script to the command
# line of the terminal.

################################################################################

# General settings
mem=32G
tmpspace=10G
nrthreads=8

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04/

################################################################################
# Mapping related settings

# Mapping parameters
genome=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/STAR-indexed-49
#gtffile=/Volumes/fastq_m.wehrens/Mapping/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/filteredcustom_Homo_sapiens.GRCh38.93.gtf
  # only protein_coding, lincRNA, antisense (like cell ranger)
  
################################################################################
# Data directory
  
parentdir=/hpc/hub_oudenaarden/mwehrens/fastq/Sjoerd
#all_fastqfiles=(head_clone-12-no-dox-rep1_AAC3KKWM5_S19_L001)
all_fastqfiles=(clone-11-dox-rep1_AAC3KKWM5_S10 clone-11-dox-rep2_AAC3KKWM5_S11 clone-11-dox-rep3_AAC3KKWM5_S12 clone-11-no-dox-rep1_AAC3KKWM5_S13 clone-11-no-dox-rep2_AAC3KKWM5_S14 clone-11-no-dox-rep3_AAC3KKWM5_S15 clone-12-dox-rep1_AAC3KKWM5_S16 clone-12-dox-rep2_AAC3KKWM5_S17 clone-12-dox-rep3_AAC3KKWM5_S18 clone-12-no-dox-rep1_AAC3KKWM5_S19 clone-12-no-dox-rep2_AAC3KKWM5_S20 clone-12-no-dox-rep3_AAC3KKWM5_S21 clone-2-dox-rep1_AAC3KKWM5_S1 clone-2-dox-rep2_AAC3KKWM5_S2 clone-2-dox-rep3_AAC3KKWM5_S3 clone-3-dox-rep1_AAC3KKWM5_S4 clone-3-dox-rep2_AAC3KKWM5_S5 clone-3-dox-rep3_AAC3KKWM5_S6 clone-5-no-dox-rep1_AAC3KKWM5_S7 clone-5-no-dox-rep2_AAC3KKWM5_S8 clone-5-no-dox-rep3_AAC3KKWM5_S9 cone-5-dox-rep1_AAC3KKWM5_S22 cone-5-dox-rep2_AAC3KKWM5_S23 cone-5-dox-rep3_AAC3KKWM5_S24)

################################################################################  
# Start the jobs

for myfastqfile in "${all_fastqfiles[@]}"
do
  
  echo "Starting jobs for $myfastqfile.."
  
  # nrthreads will be set to 1 here
  # also increased tmpspace and mem will not be requested here
  job1=$(sbatch  --parsable --output=slurm.%x.%j.out --job-name=p1_${myfastqfile} -c 1 --time=14-00:00:00  \
    --export=ALL,myfastqfile="${myfastqfile}",parentdir="${parentdir}" ${script_path}/bulk_seq_mapping/bulk_seq_mapping_part1.sh)
  
  job2=$(sbatch --dependency=afterany:${job1} --parsable --output=slurm.%x.%j.out --job-name=p2_${myfastqfile} -c ${nrthreads} --gres=tmpspace:${tmpspace} --time=14-00:00:00  --mem=${mem} \
      --export=ALL,nrthreads="${nrthreads}",genome="${genome}",myfastqfile="${myfastqfile}",gtffile="${gtffile}",parentdir="${parentdir}" ${script_path}/bulk_seq_mapping/bulk_seq_mapping_part2.sh)
        
  echo "For $myfastqfile, started job1: $job1; job2: $job2;"  
  
done










