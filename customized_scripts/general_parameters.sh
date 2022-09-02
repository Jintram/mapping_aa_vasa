#!/bin/bash

echo "Settings (server) general parameters"

# PARAMETERS FOR RUNNING A JOB ON THE HPC SERVER
p2s=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04
path2scripts=$p2s

#p2trimgalore=/hpc/hub_oudenaarden/aalemany/bin/TrimGalore-0.4.3
p2trimgalore=/hpc/hub_oudenaarden/mwehrens/bin/TrimGalore-0.6.6
path2trimgalore=$p2trimgalore

#p2cutadapt=/hpc/hub_oudenaarden/aalemany/bin
p2cutadapt=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin
path2cutadapt=$p2cutadapt
  # run "python $p2cutadapt/cutadapt --version" to check version

#p2bwa=/hpc/hub_oudenaarden/bin/software/bwa-0.7.10
p2bwa=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin
#p2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1/
p2samtools=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin

#p2star=/hpc/hub_oudenaarden/aalemany/bin/STAR-2.7.3a/bin/Linux_x86_64/
p2star=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin

pythonbin=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin/python

#p2bedtools=/hpc/hub_oudenaarden/aalemany/bin/bedtools2/bin/
p2bedtools=/hpc/hub_oudenaarden/mwehrens/bin/bedtools/
p2b=$p2bedtools





