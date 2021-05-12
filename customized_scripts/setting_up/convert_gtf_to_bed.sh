#!/bin/bash

################################################################################

# This script converts a gtf file to Anna's custom bed format

# The custom format Anna is using has the following format:
# example line:
# 1	11869	12227	+	ENSG00000223972_DDX11L1_TranscribedUnprocessedPseudogene_exon	2540	11869	14409
# The columns are:
# chromosome_nr start_pos end_pos strand ENSEMBLID_GENENAME_GENEBIOTYPE_EXONORINTRON GENE_LENGTH GENE_START GENE_END
# Note that genes are grouped together

# Execute this script, e.g. using:
# p2s=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/setting_up/
# ${p2s}/convert_gtf_to_bed.sh /Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/downloaded_ensembl/81/head_Homo_sapiens.GRCh38.81.gtf customGenestRNA_
#
# For testing purposes, execute line-by-line, using e.g.
# gtffile=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/downloaded_ensembl/81/Homo_sapiens.GRCh38.81.gtf 
# output=homeMade_IntronsExons__Homo_sapiens.GRCh38.81_

#p2scripts=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts
p2scripts=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04

if [ "$gtffile" = "" ]; then
    gtffile=$1
    output=$2
fi

python ${p2scripts}/setting_up/get_IntronsExons_fromGTF.py ${gtffile} ${output}

# after this, there's three files, one for introns, one for exons, and one for genes.
# we now still need to combine this information in such a way that the structure of the rows
# is
# ---
# exon
# intron
# exon
# intron
# exon
# etc
# ---
# And also to add the fields, GENE_LENGTH GENE_START GENE_END
# I'm guessing this should be possible using a few commands.

# Create an additional field indicating exon or intron status, and merge into _ExonsIntrons file
awk -F, '{$(NF+1)="exon";}1' OFS='\t' ${output}_mw_exons.bed > ${output}_ExonsIntrons.bed
awk -F, '{$(NF+1)="intron";}1' OFS='\t' ${output}_mw_introns.bed >> ${output}_ExonsIntrons.bed
#cat ${output}_exons.bed ${output}_introns.bed > ${output}_ExonsIntrons.bed # Paste together the intron and exon files

# Calculate gene length
awk 'BEGIN {OFS="\t"} { $(NF+1) = $3 - $2 } 1' ${output}_genes.bed > ${output}_genesL.bed

# Join exon/intron file with gene file
# echo "join requires sorted input"
sort -k5 ${output}_ExonsIntrons.bed > ${output}_ExonsIntrons_sorted.bed
sort -k5 ${output}_genesL.bed > ${output}_genesL_sorted.bed
join -t $'\t' -j 5 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.6,2.2,2.3 ${output}_ExonsIntrons_sorted.bed ${output}_genesL_sorted.bed > ${output}_joined.bed
  # DEBUG: check sorting of ENSG00000215790_SLC35E2A_ProteinCoding
  
# Sort the resulting file
sed 's/	/_/5' ${output}_joined.bed > ${output}_joined2.bed
sort -k1,1 -k2,2n <  ${output}_joined2.bed > ${output}_Final.bed 

# Remove intermediate files
rm ${output}_genesL.bed
rm ${output}_genesL_sorted.bed
rm ${output}_exons.bed
rm ${output}_introns.bed
rm ${output}_ExonsIntrons.bed
rm ${output}_ExonsIntrons_sorted.bed
rm ${output}_joined2.bed
rm ${output}_joined.bed

#${output}_mw_exons.bed
#${output}_mw_introns.bed
#${output}_genes.bed

# Let's check if consistent with Anna's file
#annasversion=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/Homosapines_ensemble99.homemade_IntronExonTrna.bed
#sort -k1,1 -k2,2n $annasversion > anna_sorted.bed 
#diff ${output}_Final.bed anna_sorted.bed > diff.txt
#diff ${output}_Final.bed $annasversion > diff.txt



################################################################################
# Some older debugging stuff

#rm XXXSTUFFXXX
#output=homeMade_IntronsExons__Homo_sapiens.GRCh38.99_
#annasversion=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/Homosapines_ensemble99.homemade_IntronExonTrna.bed
#diff ${output}_Final.bed $annasversion > diff.txt

# There seem to be some differences, perhaps she used?
#wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf.gz

################################################################################

# Old code; I started setting something up myself, but Anna quickly sent me her
# version so this was not needed any more.
  ## Add the transcript_id field if it's not present, and continue with conversion
  #awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' head_Homo_sapiens.GRCh38.81.gtf | gtf2bed - > head_Homo_sapiens.GRCh38.81.bed
  ##gtf2bed < head_Homo_sapiens.GRCh38.81.gtf
  #
  #python convert_gtf_to_bed_implement.py
  
################################################################################