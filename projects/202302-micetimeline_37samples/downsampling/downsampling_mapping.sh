#!/bin/bash

# This script will perform some mapping on downsampled data, to check whether
# additional sequencing might yield additional results.

# We'll start working with the fastq files for which already pT reads
# were selected, and R1 and R2 were already combined into 1 file, such that
# BCs etc are known. (--> Skip step 1 of mapping.)

# My assumption is that the reads are in random order, such that downsampling
# can simply be done with the head command.

# The total amount of reads in the file are:
# 30205212/4 = 7551303 (wc -l command)

# Let's simply create 10%, 20%, 30% etc. files
# 10% ± equals 3020520 (deliberately a multiple of 4)

# To start interactive node:
# srun --nodes=1 --time=02:00:00 --pty bash -i

# ==============================================================================
# Perform downsampling

ds_values=("4000" "3020520" "6041040" "9061560" "12082080" "15102600" "18123120" "21143640" "24164160" "27184680" "30205200")

for i in ${ds_values[@]}; do
  echo "Creating file with $i lines.."
  zcat Martijn-Wehrens-sample1_AACGFYKM5_S2_L001_pT_R2_cbc.fastq.gz | head -n $i > head${i}__Martijn-Wehrens-sample1_AACGFYKM5_S2_L001_pT_R2_cbc.fastq
  
  gzip head${i}__Martijn-Wehrens-sample1_AACGFYKM5_S2_L001_pT_R2_cbc.fastq
done
  
# ==============================================================================  
# Now start the jobs

# Repeat for convenience
ds_values=("4000" "3020520" "6041040" "9061560" "12082080" "15102600" "18123120" "21143640" "24164160" "27184680" "30205200")

# Paths
script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
general_parameter_filepath=${script_path}/general_parameters.sh
run_parameter_filepath=/hpc/hub_oudenaarden/mwehrens/fastq/202302_micetimeline_37samples/run_parameters_mice37_genome-n-ERCC_downsampling.sh
  # (Step 1 will be skipped)

# Start jobs
for i in ${ds_values[@]}; do
  
  # i="3020520"
  # i=4000
  
  lib=head${i}__Martijn-Wehrens-sample1_AACGFYKM5_S2_L001 # Fastq files name, without _Rx.fastq.gz
  tmpspace=10G
  mymem=50G
  loadtmpwithoutdir=1
  
  echo "Starting lib $lib"
  
  sbatch --job-name=s1dsmice --time=72:00:00   --gres=tmpspace:${tmpspace} --mem=${mymem} --export=ALL,loadtmpwithoutdir="${loadtmpwithoutdir}",general_parameter_filepath="${general_parameter_filepath}",run_parameter_filepath="${run_parameter_filepath}",lib="$lib" ${script_path}/L_TS_submit_vasaplate_map_LOCAL.sh

done

# ==============================================================================  
# Now create count tables..


mappingfolder=/hpc/hub_oudenaarden/mwehrens/fastq/202302_micetimeline_37samples/mapping_downsampletest/

cd $mappingfolder

mem=10G
tmpspace=10G

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04/
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107.ERCC/GRCm39.dna.primary_plus_ERCC.gtf
  # only protein_coding, lincRNA, antisense (like cell ranger)
  
trimUMI=0
countmulti=no

ds_values=("4000" "3020520" "6041040" "9061560" "12082080" "15102600" "18123120" "21143640" "24164160" "27184680" "30205200")

for i in ${ds_values[@]}
do
  #sample=$s1
  sample=head${i}__Martijn-Wehrens-sample1_AACGFYKM5_S2
  for lib in ${sample}_L001_pT 
  do
    echo $lib
    
    bamfile=${mappingfolder}/${lib}.nonRibo_E99_Aligned.out
    # bamfile=/Volumes/fastq_m.wehrens/Mapping/WANG4/mapping/${lib}_nc.nonRibo_E99_Aligned.out
    
    nrthreads=8
    job1=$(sbatch  --parsable --output=slurm-${lib}-%x.%j.out --job-name=cntUT_${lib} -c ${nrthreads} --gres=tmpspace:${tmpspace} --time=14-00:00:00  --mem=${mem} \
      --export=ALL,gtffile="${gtffile}",nrthreads="${nrthreads}",trimUMI="${trimUMI}",steptorun="1",bamfile="${bamfile}",countmulti="${countmulti}" ${script_path}/countTables_umiTools.sh p.N1.plate.97474.part.1_cat)
    
    nrthreads=1
    job2=$(sbatch --dependency=afterany:${job1} --parsable --output=slurm-${lib}-%x.%j.out --job-name=cntUT_${lib} -c ${nrthreads} --gres=tmpspace:${tmpspace} --time=14-00:00:00  --mem=${mem} \
        --export=ALL,gtffile="${gtffile}",nrthreads="${nrthreads}",trimUMI="${trimUMI}",steptorun="2",bamfile="${bamfile}",countmulti="${countmulti}" ${script_path}/countTables_umiTools.sh p.N1.plate.97474.part.1_cat)
    
    nrthreads=1
    job3=$(sbatch --dependency=afterany:${job2} --parsable --output=slurm-${lib}-%x.%j.out --job-name=cntUT_${lib} -c ${nrthreads} --gres=tmpspace:${tmpspace} --time=14-00:00:00  --mem=${mem} \
        --export=ALL,gtffile="${gtffile}",nrthreads="${nrthreads}",trimUMI="${trimUMI}",steptorun="3",bamfile="${bamfile}",countmulti="${countmulti}" ${script_path}/countTables_umiTools.sh p.N1.plate.97474.part.1_cat)
    
      
    echo "started job1: $job1; job2: $job2; job3: $job3"  
    
  done
done








