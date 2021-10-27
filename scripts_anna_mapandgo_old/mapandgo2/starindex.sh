#!/bin/bash

#genome=hpc/hub_oudenaarden/group_references/ensembl/93/homo_sapiens
#genome=/hpc/hub_oudenaarden/group_references/ensembl/94/c_elegans
genome=/hpc/hub_oudenaarden/group_references/ensembl/93/danio_rerio
genomeDir=${genome}/star_index_75
genomeFastq=${genome}/primary_assembly_ERCC92.fa
genomeGTF=${genome}/annotations.gtf

path2star=/hpc/hub_oudenaarden/avo/nascent/STAR-2.5.3a/bin/Linux_x86_64

${path2star}/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFiles ${genomeFastq} --sjdbGTFfile ${genomeGTF} --sjdbOverhang 75


