#!/bin/bash

# For bulk mapping without UMIs or BCs, the protocol is very simple:
# 
# 1. Perform trimming of the fastq files
# 2. Use STAR to map the files
# 3. Use featureCounts to convert to count table
#
# Step 2 & 3 can be performed multi-threaded, so these scripts are split in 
# two parts, (A) step 1, (B) step 2+3.

# For testing purposes:
# srun --nodes=1 --time=02:00:00 --mem=50G --pty bash -i 
# myfastqfile=head_clone-12-no-dox-rep1_AAC3KKWM5_S19_L001

# Paths
path2trimgalore=/hpc/hub_oudenaarden/mwehrens/bin/TrimGalore-0.6.6
path2cutadapt=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin

echo "Trimming; $myfastqfile"

# Let's do the trimming:
${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt --paired ${myfastqfile}_R1.fastq.gz ${myfastqfile}_R2.fastq.gz

# This result in the files 
# ${myfastqfile}_R1_val_1.fq.gz
# ${myfastqfile}_R2_val_2.fq.gz

echo "Trimming done ..."






