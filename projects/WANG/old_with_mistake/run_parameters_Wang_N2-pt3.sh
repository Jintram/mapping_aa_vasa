#!/bin/bash

echo "Setting (server) run parameters"

filetorunbefore=wang_n2_continue.sh

# format for step: step="(1)(3)" etc or step="default", latter runs all default
# note that at the server, temporary files might be deleted ..
step="(5)(6)"

#nocleanup=1
#returntempfiles=1 # on HPC, copy tempfiles from $TMPDIR to $outdir

outdir=/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/mapping_may25

# read_length - 1 (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136)
n=49

# human genome
#riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/unique_rRNA_human.fa
#genome=/hpc/hub_oudenaarden/group_references/ensembl/99/homo_sapiens/star_v273a_NOMASK_NOERCC_index_$n
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Homosapines_ensemble99.homemade_IntronExonTrna.bed
riboref=/hpc/hub_oudenaarden/mwehrens/ref/rRNA/unique_rRNA_human.fa
genome=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/STAR-indexed-49
refBED=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/homeMadeBed/homeMade_IntronsExons__Homo_sapiens.GRCh38.93__Final.bed

# mouse genome
#riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
#genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

# choose protocol from [celseq1, celseq2, vasaplate]
protocol="takara"

# Whether we should look for targeted reads
TS=0

stranded=y







