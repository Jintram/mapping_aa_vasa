#!/bin/bash

echo "Setting (server) run parameters"

# format for step: step="(1)(3)" etc or step="default", latter runs all default
step="default"

outdir=/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping.93.may25

# read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136)
n=74 

# human genome
#riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/unique_rRNA_human.fa
#genome=/hpc/hub_oudenaarden/group_references/ensembl/99/homo_sapiens/star_v273a_NOMASK_NOERCC_index_$n
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Homosapines_ensemble99.homemade_IntronExonTrna.bed
riboref=/hpc/hub_oudenaarden/mwehrens/ref/rRNA/unique_rRNA_human.fa
  # I assume /pub/release-93/fasta/homo_sapiens/ncrna , filter for biotype (but this file was obtained from anna)
genome=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/STAR-indexed-74
  # Note: this can be downloaded from 
  # - ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  
  # - ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf
refBED=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/homeMadeBed/homeMade_IntronsExons__Homo_sapiens.GRCh38.93__Final.bed


# mouse genome
#riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
#genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

# choose protocol from [celseq1, celseq2, vasaplate]
protocol="celseq2"

# Whether we should look for targeted reads
TS=0

stranded=y




