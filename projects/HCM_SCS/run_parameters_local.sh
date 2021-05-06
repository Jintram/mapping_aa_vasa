#!/bin/bash

# 
echo "Setting (local) run parameters"

step="2"

outdir=/Volumes/workdrive_m.wehrens_hubrecht/Fastq__raw_files/HCM_SCS/testrun/mapping/

# read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136)
n=74 # note this has to do with the pre-processing of the reference genome

# human genome
riboref=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/unique_rRNA_human.fa
genome=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/star_v273a_NOMASK_NOERCC_index_74    
refBED=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/Homosapines_ensemble99.homemade_IntronExonTrna.bed

# mouse genome
#riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
#genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

# choose protocol from [celseq1, celseq2, vasaplate]
protocol="celseq2"

# Whether we should look for targeted reads
TS=0
