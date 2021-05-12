#!/bin/bash

# Another way to get dir of scripts
#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

### check input parameters
#if [[ $# -ne 2 ]] && [[ -z "$general_parameter_filepath" || -z "$run_parameter_filepath" || -z "$lib"]]
if [[ ($# -ne 3) && (-z "$general_parameter_filepath" || -z "$run_parameter_filepath" || -z "$lib") ]]
then
    echo "Please, give location of files with general and run parameters or set general_parameter_filepath and run_parameter_filepath."
    echo "1) general_parameter_filepath"
    echo "2) run_parameter_filepath"
    echo "3) library prefix (filename of fastQ without _R1.fastq.gz)"
    exit
fi

# Set paths
# general
if [ -z "$general_parameter_filepath" ]
then
  general_parameter_filepath=$1
fi
# run-specific
if [ -z "$run_parameter_filepath" ]
then
  run_parameter_filepath=$2
fi
# lib
if [ -z "$lib" ]
then
  lib=$3
fi

echo "Running script on ${lib} library"
echo "Running step $step"

# Load parameters from paths
# general
if [ ! -f "$general_parameter_filepath" ]
then
    echo "$general_parameter_filepath not found"
    exit
fi
echo "sourcing $general_parameter_filepath"
source $general_parameter_filepath
# run-specific
if [ ! -f "$run_parameter_filepath" ]
then
    echo "$run_parameter_filepath not found"
    exit
fi
echo "sourcing $run_parameter_filepath"
source $run_parameter_filepath

if [ "$step" = "" ] 
then
  echo "Step was not set, setting to default (full run)"
  step="default"
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
if [[ "${step}" == *"(1)"* || "${step}" == "default" ]]; then
    sh ${p2s}/L_TS_extractBC.sh $general_parameter_filepath $run_parameter_filepath vasaplate $lib
fi

if [[ "${step}" == *"(2)"* || "${step}" == "default" ]]; then
    ### trim files
    # (no local version needed)
    # For polyT:
    ${p2s}/TS_trim.sh $general_parameter_filepath $run_parameter_filepath ${lib}_pT_R2
    # For non-classified:
    ${p2s}/TS_trim.sh $general_parameter_filepath $run_parameter_filepath ${lib}_nc_R2 
    # For TS
    if [ $TS = '1' ]; then
      ${p2s}/TS_trim_paired.sh $general_parameter_filepath $run_parameter_filepath ${lib}_TS 
    fi
fi

if [[ "${step}" == *"(3)"* || "${step}" == "default" ]]; then

      echo "Mapping to ribosomal RNA. (More relevant for VASA.)"
  
      ### ribo-map
      # this is done "in silico" remove ribosomal reads (should be filtered by wet lab protocol already, but not 100%)
      # not that ribo-bwamem script also removes the reads using riboread-selection.py after mapping
      # 
      # map poly-T reads to ribo
      ${p2s}/L_TS_ribo-bwamem.sh $general_parameter_filepath $run_parameter_filepath ${lib}_pT_R2_cbc_trimmed_HATCG.fq.gz ${lib}_pT_cbc_trimmed_HATCG 
      # map unclassified reads to ribo
      ${p2s}/L_TS_ribo-bwamem.sh $general_parameter_filepath $run_parameter_filepath ${lib}_nc_R2_cbc_trimmed_HATCG.fq.gz ${lib}_nc_cbc_trimmed_HATCG 
      # map targeted reads to ribo
      if [ $TS = '1' ]; then
        ${p2s}/L_TS_ribo-bwamem_paired.sh ${lib}_TS_R1_cbc_val_1_HATCG.fq.gz ${lib}_TS_R2_cbc_val_2_HATCG.fq.gz ${lib}_TS_cbc_val_HATCG $p2bwa $p2samtools y $p2s
      fi
    
fi      
    
if [[ "${step}" == *"(4)"* || "${step}" == "default" ]]; then
    ### map to genome 
    # For the poly-T data (single)
    ${p2s}/L_TS_map_star.sh $general_parameter_filepath $run_parameter_filepath ${lib}_pT_cbc_trimmed_HATCG.nonRibo.fastq ${lib}_pT.nonRibo_E99_
    # Then for the unclassified data (single)
    ${p2s}/L_TS_map_star.sh $general_parameter_filepath $run_parameter_filepath ${lib}_nc_cbc_trimmed_HATCG.nonRibo.fastq ${lib}_nc.nonRibo_E99_
    # For the paired TS data:
    if [ $TS = '1' ]; then      
      ${p2s}/L_TS_map_star_paired.sh $general_parameter_filepath $run_parameter_filepath ${lib}_TS_cbc_val_HATCG_R1.nonRibo.fastq ${lib}_TS_cbc_val_HATCG_R2.nonRibo.fastq ${lib}_TS.nonRibo_E99_
        # note: locally, .gz is removed from input file names
    fi
fi

if [[ "${step}" == *"(5)"* || "${step}" == "default" ]]; then
    ### map locations to genes (accounting for ambiguities)
      # note that there can be ambiguities in where stuff needs to be mapped
      # mostly due to overlapping annotation for whatever reason or
      # multi-mappers
      #
      # this script uses awk commands to apply a set of pre-determined rules for 
      # dealing with ambiguities in the assignment of locations to ref transcripts
    #
    # For pT subset
    ${p2s}/L_TS_deal_with_singleandmultimappers_paired_stranded.sh  $general_parameter_filepath $run_parameter_filepath ${lib}_pT.nonRibo_E99_Aligned.out.bam n
    # For nc subset
    ${p2s}/L_TS_deal_with_singleandmultimappers_paired_stranded.sh  $general_parameter_filepath $run_parameter_filepath ${lib}_nc.nonRibo_E99_Aligned.out.bam n
    # For TS subset
    if [ $TS = '1' ]; then
      ${p2s}/L_TS_deal_with_singleandmultimappers_paired_stranded.sh  $general_parameter_filepath $run_parameter_filepath ${lib}_TS.nonRibo_E99_Aligned.out.bam y
    fi
      
fi

if [[ "${step}" == *"(6)"* || "${step}" == "default" ]]; then
      
    ### count table    
    # For pT
    $pythonbin ${p2s}/TS_countTables_final.py ${outdir}/${lib}_pT.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${outdir}/${lib}_pT.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${outdir}/${lib}_pT vasa
    # For nc
    $pythonbin ${p2s}/TS_countTables_final.py ${outdir}/${lib}_nc.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${outdir}/${lib}_nc.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${outdir}/${lib}_nc vasa
    # For TS
    if [ $TS = '1' ]; then
      $pythonbin ${p2s}/TS_countTables_final.py ${outdir}/${lib}_TS.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${outdir}/${lib}_TS.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${outdir}/${lib}_TS vasa 1
    fi

fi

if [[ "${step}" == *"(7)"* ]]; then  

    ### TS specific post-processing, to get a database of mapped reads and their sequences
    # TS_countTables_final should output "decision" files also, which contains final mapping decisions for multimappers
    # (single mappers also included for easy processing)
    # First merge these files, then sort them
    if [ $TS = '1' ]; then
      cat ${lib}_singlemapper_decisions.tsv > ${lib}_merged_decisions.tsv
      cat ${lib}_multimapper_decisions.tsv >> ${lib}_merged_decisions.tsv
      sort ${lib}_merged_decisions.tsv > ${lib}_merged_decisions_sorted.tsv
      ${p2samtools}/samtools sort -n -O bam -T temp_${lib} -o ${lib}_TS.nonRibo_E99_Aligned.nsorted.out.bam ${lib}_TS.nonRibo_E99_Aligned.out.bam
      $pythonbin ${p2s}/TS_mw_sequence_db.py ${lib}
    fi
    
    # for debugging purposes, bam->sam
    # samtools view -h -o ${lib}_TS.nonRibo_E99_Aligned.nsorted.out.sam ${lib}_TS.nonRibo_E99_Aligned.nsorted.out.bam

    # done
fi








