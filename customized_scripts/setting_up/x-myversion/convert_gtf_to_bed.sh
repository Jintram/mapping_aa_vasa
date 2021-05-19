#!/bin/bash


# This script converts a gtf file to Anna's custom bed format

# The custom format Anna is using has the following format:
# example line:
# 1	11869	12227	+	ENSG00000223972_DDX11L1_TranscribedUnprocessedPseudogene_exon	2540	11869	14409
# The columns are:
# chromosome_nr start_pos end_pos strand ENSEMBLID_GENENAME_GENEBIOTYPE_EXONORINTRON GENE_LENGTH GENE_START GENE_END
# Note that genes are grouped together

# Add the transcript_id field if it's not present, and continue with conversion
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' head_Homo_sapiens.GRCh38.81.gtf | gtf2bed - > head_Homo_sapiens.GRCh38.81.bed

#gtf2bed < head_Homo_sapiens.GRCh38.81.gtf

python convert_gtf_to_bed_implement.py