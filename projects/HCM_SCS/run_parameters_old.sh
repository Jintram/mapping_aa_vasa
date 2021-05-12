#!/bin/bash

echo "Setting (server) run parameters"

# format for step: step="(1)(3)" etc or step="default", latter runs all default
step="(5)"

outdir=/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping

# read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136)
n=74 

# human genome
riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/unique_rRNA_human.fa
genome=/hpc/hub_oudenaarden/group_references/ensembl/99/homo_sapiens/star_v273a_NOMASK_NOERCC_index_$n
refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Homosapines_ensemble99.homemade_IntronExonTrna.bed

# mouse genome
#riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
#genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

# choose protocol from [celseq1, celseq2, vasaplate]
protocol="celseq2"

# Whether we should look for targeted reads
TS=0

stranded=y
paired=0