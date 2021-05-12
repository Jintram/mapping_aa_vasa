#!/bin/bash

# submit by:
# sbatch --job-name=bwaIndx --time=12:00:00  --mem=100G ${script_path}/setting_up/create-rRNA-index.sh

bwa index /hpc/hub_oudenaarden/mwehrens/ref/rRNA/unique_rRNA_human.fa