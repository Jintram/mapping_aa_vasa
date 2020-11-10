#!/bin/bash

# LOCAL VERSION

### check input parameters
if [ $# -ne 4 ] 
then
    echo "Please, give:"
    echo "1) library name (prefix of the fastq files)"
    echo "2) genome: MOUSE /  HUMAN"
    echo "3) read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136)"
    echo "4) prefix for output files"
    exit
fi

lib=$1
ref=$2
n=$3
out=$4

### input paths (to modify by user)

#p2s=/home/hub_oudenaarden/mwehrens/scripts_AAMW
#p2trimgalore=/hpc/hub_oudenaarden/aalemany/bin/TrimGalore-0.4.3
#p2cutadapt=/hpc/hub_oudenaarden/aalemany/bin
#p2bwa=/hpc/hub_oudenaarden/bin/software/bwa-0.7.10
#p2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1/
#p2star=/hpc/hub_oudenaarden/aalemany/bin/STAR-2.7.3a/bin/Linux_x86_64/
#p2bedtools=/hpc/hub_oudenaarden/aalemany/bin/bedtools2/bin/

p2s=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts
p2trimgalore=/Users/m.wehrens/Software_custom/TrimGalore-0.6.6
p2cutadapt=/Users/m.wehrens/Library/Python/3.7/bin
p2bwa=/Users/m.wehrens/Software_custom/bwa
p2samtools=/Users/m.wehrens/Software_custom/samtools-1.2
p2star=/usr/local/bin
pythonbin=/Users/m.wehrens/anaconda3/bin/python

### check existence of input fastq files
if [ ! -f ${lib}_R1.fastq.gz ]
then
    echo "${lib}_R1.fastq.gz not found"
    exit
fi

if [ ! -f ${lib}_R2.fastq.gz ]
then
    echo "${lib}_R2.fastq.gz not found"
    exit
fi

### check python version (we want version 3)
v=$($pythonbin -c 'import sys; print(".".join(map(str, sys.version_info[:3])))' | awk -F "." '{print $1}')
if [ $v -ne "3" ]
then
    echo "python needs to be 3"
    exit
fi

### set references
if [[ $ref == "MOUSE" ]]
then
    riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
    genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
    if [ ! -d $genome ]
    then
        echo "genome not found"
        exit
    fi
    refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
elif [[ $ref == "HUMAN" ]]
then
    #riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/unique_rRNA_human.fa
    #genome=/hpc/hub_oudenaarden/group_references/ensembl/99/homo_sapiens/star_v273a_NOMASK_NOERCC_index_$n
    riboref=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/unique_rRNA_human.fa
    genome=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/star_v273a_NOMASK_NOERCC_index_74    
    if [ ! -d $genome ]
    then
        echo "genome not found"
        exit
    fi
    refBED=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/Homosapines_ensemble99.homemade_IntronExonTrna.bed
    #refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Homosapines_ensemble99.homemade_IntronExonTrna.bed
fi

### extract cell barcodes (this will split the read files into 3 batches; poly-T, targeted and unclassified.)
${p2s}/L_TS_extractBC.sh ${lib} vasaplate ${p2s}
  
### trim files
# (no local version needed)
# For TS
${p2s}/TS_trim_paired.sh ${lib}_TS ${p2trimgalore} ${p2cutadapt}
# For polyT:
${p2s}/TS_trim.sh ${lib}_pT_R2 ${p2trimgalore} ${p2cutadapt}
# For non-classified:
${p2s}/TS_trim.sh ${lib}_nc_R2 ${p2trimgalore} ${p2cutadapt}

### ribo-map
# map poly-T reads to ribo
${p2s}/L_TS_ribo-bwamem.sh $riboref ${lib}_pT_R2_cbc_trimmed_HATCG.fq.gz ${lib}_pT_cbc_trimmed_HATCG $p2bwa $p2samtools y $p2s
# map targeted reads to ribo
${p2s}/L_TS_ribo-bwamem_paired.sh $riboref ${lib}_TS_R1_cbc_val_1_HATCG.fq.gz ${lib}_TS_R2_cbc_val_2_HATCG.fq.gz ${lib}_TS_cbc_val_HATCG $p2bwa $p2samtools y $p2s
# map unclassified reads to ribo
${p2s}/L_TS_ribo-bwamem.sh $riboref ${lib}_nc_R2_cbc_trimmed_HATCG.fq.gz ${lib}_nc_cbc_trimmed_HATCG $p2bwa $p2samtools y $p2s
  # this is done "in silico" remove ribosomal reads (should be filtered by wet lab protocol already, but not 100%)
  # not that ribo-bwamem script also removes the reads using riboread-selection.py after mapping

### map to genome 
# First for the paired TS data:
${p2s}/L_TS_map_star_paired.sh ${p2star} ${p2samtools} ${genome} ${lib}_TS_cbc_val_HATCG_R1.nonRibo.fastq ${lib}_TS_cbc_val_HATCG_R2.nonRibo.fastq ${lib}_TS.nonRibo_E99_
  # note: locally, .gz is removed from input file names
# Then for the poly-T data (single)
${p2s}/L_TS_map_star.sh ${p2star} ${p2samtools} ${genome} ${lib}_pT_cbc_trimmed_HATCG.nonRibo.fastq ${lib}_pT.nonRibo_E99_
# Then for the unclassified data (single)
${p2s}/L_TS_map_star.sh ${p2star} ${p2samtools} ${genome} ${lib}_nc_cbc_trimmed_HATCG.nonRibo.fastq ${lib}_nc.nonRibo_E99_


### map locations to genes (accounting for ambiguities)
  # note that there can be ambiguities in where stuff needs to be mapped
  # mostly due to overlapping annotation for whatever reason or
  # multi-mappers
  
# For TS subset
${p2s}/L_TS_deal_with_singlemappers_paired.sh ${lib}_TS.nonRibo_E99_Aligned.out.bam ${refBED} y
# For pT subset
${p2s}/L_TS_deal_with_singlemappers_paired.sh ${lib}_pT.nonRibo_E99_Aligned.out.bam ${refBED} y
# For nc subset
${p2s}/L_TS_deal_with_singlemappers_paired.sh ${lib}_nc.nonRibo_E99_Aligned.out.bam ${refBED} y
  # this script uses awk commands to apply a set of pre-determined rules for 
  # dealing with ambiguities in the assignment of locations to ref transcripts

# For TS subset
${p2s}/L_TS_deal_with_multimappers_paired.sh ${lib}_TS.nonRibo_E99_Aligned.out.bam ${refBED} y
# For pT subset
${p2s}/L_TS_deal_with_multimappers_paired.sh ${lib}_pT.nonRibo_E99_Aligned.out.bam ${refBED} y
# For nc subset
${p2s}/L_TS_deal_with_multimappers_paired.sh ${lib}_nc.nonRibo_E99_Aligned.out.bam ${refBED} y


### count table
# For TS
$pythonbin ${p2s}/TS_countTables_final.py ${lib}_TS.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${lib}_TS.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${lib}_TS vasa
# For pT
$pythonbin ${p2s}/TS_countTables_final.py ${lib}_pT.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${lib}_pT.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${lib}_pT vasa
# For nc
$pythonbin ${p2s}/TS_countTables_final.py ${lib}_nc.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${lib}_nc.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${lib}_nc vasa

### TS specific post-processing, to get a database of mapped reads and their sequences
# TS_countTables_final should output "decision" files also, which contains final mapping decisions for multimappers
# (single mappers also included for easy processing)
# First merge these files, then sort them
cat ${lib}_singlemapper_decisions.tsv > ${lib}_merged_decisions.tsv
cat ${lib}_multimapper_decisions.tsv >> ${lib}_merged_decisions.tsv
sort ${lib}_merged_decisions.tsv > ${lib}_merged_decisions_sorted.tsv
${p2samtools}/samtools sort -n -O bam -T temp_${lib} -o ${lib}_TS.nonRibo_E99_Aligned.nsorted.out.bam ${lib}_TS.nonRibo_E99_Aligned.out.bam
$pythonbin ${p2s}/TS_mw_sequence_db.py ${lib}

# for debugging purposes, bam->sam
# samtools view -h -o ${lib}_TS.nonRibo_E99_Aligned.nsorted.out.sam ${lib}_TS.nonRibo_E99_Aligned.nsorted.out.bam

# done









