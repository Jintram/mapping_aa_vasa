#!/bin/bash

echo "Setting (server) run parameters"

# format for step: step="(1)(3)" etc or step="default", latter runs all default
step="(1)(2)(3)(4)"

#outdir=/hpc/hub_oudenaarden/mwehrens/fastq/202209_micetimeline/mapping
outdir=/hpc/hub_oudenaarden/mwehrens/fastq/202308_HFTFs_mapping

# read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136)
n=100 
  # note this needs to be (read length)-1

# human genome
##riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/unique_rRNA_human.fa
##genome=/hpc/hub_oudenaarden/group_references/ensembl/99/homo_sapiens/star_v273a_NOMASK_NOERCC_index_$n
##refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Homosapines_ensemble99.homemade_IntronExonTrna.bed
riboref=/hpc/hub_oudenaarden/mwehrens/ref/rRNA/unique_rRNA_human.fa
genome=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.110.ERCC/STAR-indexed-100
#refBED=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.81/homeMadeBed/homeMade_IntronsExons__Homo_sapiens.GRCh38.81__Final.bed

# mouse genome
##riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
##genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
##refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
#riboref=/hpc/hub_oudenaarden/mwehrens/ref/rRNA/rRNA_mouse_Rn45S.fa
#genome=/hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107.ERCC/STAR-indexed-$n
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
  # Only necessary for custom annotation

# choose protocol from [celseq1, celseq2, vasaplate]
protocol="celseq2"

# Whether we should look for targeted reads
TS=0

stranded=y




