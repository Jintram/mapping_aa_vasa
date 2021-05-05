#!/bin/bash

# LOCAL VERSION

### check input parameters
if [ $# -ne 1 ] 
then
    echo "Please, give location of file with run parameters."
    exit
fi

# Set run-specific parameters
run_parameter_filepath=$1
if [ ! -f "$1" ]
then
    echo "$1 not found"
    exit
fi
source ./set_paths.sh

if [ "$step" = "" ] || [ "$step" = "all" ]; then
  echo "Step was not set, setting to all"
  step="all"
fi

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

### Check whether genome files are there
if [ ! -d $genome ]
then
    echo "genome not found"
    exit
fi

### extract cell barcodes (this will split the read files into 3 batches; poly-T, targeted and unclassified.)
if [ "$step" = "1" ] || [ "$step" = "all" ]; then
    ${p2s}/L_TS_extractBC.sh $run_parameter_filepath vasaplate
fi
  
if [ "$step" = "2" ] || [ "$step" = "all" ]; then
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
    ${p2s}/L_TS_deal_with_singleandmultimappers_paired_stranded.sh ${lib}_TS.nonRibo_E99_Aligned.out.bam ${refBED} y y
    # For pT subset
    ${p2s}/L_TS_deal_with_singleandmultimappers_paired_stranded.sh ${lib}_pT.nonRibo_E99_Aligned.out.bam ${refBED} y n
    # For nc subset
    ${p2s}/L_TS_deal_with_singleandmultimappers_paired_stranded.sh ${lib}_nc.nonRibo_E99_Aligned.out.bam ${refBED} y n
      # this script uses awk commands to apply a set of pre-determined rules for 
      # dealing with ambiguities in the assignment of locations to ref transcripts

    ### count table
    # For TS
    $pythonbin ${p2s}/TS_countTables_final.py ${lib}_TS.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${lib}_TS.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${lib}_TS vasa 1
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
fi








