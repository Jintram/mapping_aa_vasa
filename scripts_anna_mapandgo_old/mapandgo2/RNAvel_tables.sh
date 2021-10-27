#!/bin/bash

path2star=/hpc/hub_oudenaarden/avo/nascent/STAR-2.5.3a/bin/Linux_x86_64
#starMouseRef=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/star_index_75
if [ $# -ne 3 ]
then
    echo "Please, give (1) input bamfile; (2) genome [mouse, human, zebrafish, spiny]; (3) root for output file"
    exit
fi

bamfile=$1
out=$3

if [ $2 == 'human' ]
then
    intron=/hpc/hub_oudenaarden/group_references/ensembl/93/homo_sapiens/annotations_hs_introns_exonssubtracted.bed
    exon=/hpc/hub_oudenaarden/group_references/ensembl/93/homo_sapiens/annotations_hs_exons.bed
elif [ $2 == 'mouse' ] 
then
    intron=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/annotations_ensembl_93_mm_introns_exonsubtracted.bed
    exon=/hpc/hub_oudenaarden/group_references/ensembl/93/mus_musculus/annotations_ensembl_93_mm_exons.bed
elif [ $2 == 'zebrafish' ]
then
    intron=/hpc/hub_oudenaarden/group_references/ensembl/93/danio_rerio/zebrafish_ensemble93_introns_exonsubtracted.bed
    exon=/hpc/hub_oudenaarden/group_references/ensembl/93/danio_rerio/zebrafish_ensemble93_exons.bed
elif [ $2 == 'spiny' ]
then
    intron=/hpc/hub_oudenaarden/aalemany/henriette/ref_genome_spiny/introns.bed
    exon=/hpc/hub_oudenaarden/aalemany/henriette/ref_genome_spiny/exons.bed
fi

/hpc/hub_oudenaarden/aalemany/bin/RNAvelocity/getIntronsExons.sh ${bamfile} $intron $exon ${out}
