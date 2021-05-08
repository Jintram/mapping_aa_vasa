#!/bin/sh
#SBATCH -o ./slurm_out/slurm.%N.%j.out # STDOUT
#SBATCH -e ./slurm_out/slurm.%N.%j.err # STDERR

sbatch --time=05:15:00  --mem=10G mw_generate_indices_ranger.sh
