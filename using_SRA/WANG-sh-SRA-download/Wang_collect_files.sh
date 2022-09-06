#!/bin/bash

# To submit to HPC, use:
# sbatch --job-name=B52 -c 2 --time=48:00:00  --mem=20G Wang_collect_files.sh
#
# Or,
# td="GSM2970361_N1_LV"
# td="GSM3449620_N14"
# td="GSM2970358_N2_LV GSM2970362_N3_LV GSM2970366_N4_LV";tdn="GSM2970358_N2_LV-GSM2970362_N3_LV-GSM2970366_N4_LV"
# td="GSM2970360_N5_LV GSM2970359_N6_LA GSM2970369_N7_LA";tdn="GSM2970360_N5_LV-GSM2970359_N6_LA-GSM2970369_N7_LA"
# td="GSM2970368_N8_LA GSM2970365_N9_LA GSM2970363_N10_LA";tdn="GSM2970368_N8_LA-GSM2970365_N9_LA-GSM2970363_N10_LA"
# td="GSM2970364_N11_LA GSM2970367_N12_LA GSM3449619_N13";tdn="GSM2970364_N11_LA-GSM2970367_N12_LA-GSM3449619_N13"
# sbatch -c 2 --output=slurm-${tdn}-%x.%j.out --job-name=B52 --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${td}" Wang_collect_files.sh

# Already done:
# 
# To do:

# GSM2970361_N1_LV GSM2970358_N2_LV GSM2970362_N3_LV GSM2970366_N4_LV GSM2970360_N5_LV GSM2970359_N6_LA GSM2970369_N7_LA GSM2970368_N8_LA GSM2970365_N9_LA GSM2970363_N10_LA GSM2970364_N11_LA GSM2970367_N12_LA GSM3449619_N13 GSM3449620_N14
if [[ "$INPUT_IDENTIFIERS" == "" ]]; then
  INPUT_IDENTIFIERS=$1
fi
if [[ "$INPUT_IDENTIFIERS" == "" ]]; then
  echo "No identifier given. Exiting .."
  exit
fi

datapath='/hpc/hub_oudenaarden/mwehrens/fastq/WANG4/'
# datapath='/hpc/hub_oudenaarden/mwehrens/fastq/WANG/'
#datapath=/Volumes/workdrive_m.wehrens_hubrecht/Fastq__raw_files/Wang/

for IDENTIFIER in $INPUT_IDENTIFIERS 
#for IDENTIFIER in GSM2970361_N1_LV GSM2970358_N2_LV GSM2970362_N3_LV GSM2970366_N4_LV GSM2970360_N5_LV GSM2970359_N6_LA GSM2970369_N7_LA GSM2970368_N8_LA GSM2970365_N9_LA GSM2970363_N10_LA GSM2970364_N11_LA GSM2970367_N12_LA GSM3449619_N13 GSM3449620_N14
do

      check_FLAG=1

      cd $datapath/${IDENTIFIER}/

      filename=${IDENTIFIER}_SRR_Acc_List.txt

      echo "Collecting data for ${IDENTIFIER}.."
      # clear / open the files
      echo -n "" > ../${IDENTIFIER}_cat_R1.fastq
      echo -n "" > ../${IDENTIFIER}_cat_R2.fastq

      while read line; do
        
          # check read 1 file
          fq1file="${line}_1.fastq"
          if [ -f "${fq1file}" ]; then              
              cat ${fq1file} >> ../${IDENTIFIER}_cat_R1.fastq
          else 
              echo "$fq1file WAS MISSING! OPERATION ABORTED"
              check_FLAG=0
              exit
          fi
          
          # check read 1 file
          fq2file="${line}_2.fastq"
          if [ -f "${fq2file}" ]; then              
              cat ${fq2file} >> ../${IDENTIFIER}_cat_R2.fastq
          else 
              check_FLAG=0
              echo "$fq2file WAS MISSING! OPERATION ABORTED"
              exit
          fi
        
          echo -n "."
        
      done < $filename

      echo ""

      if [[ ${check_FLAG} == 1 ]]; then
          echo "All files found (and hopefully collected)."
      else
          echo "FILES WERE MISSING! OPERATION ABORTED"
          exit # note these lines are redundant, exit will have occured already
      fi

      echo "Zipping files"
      gzip ../${IDENTIFIER}_cat_R1.fastq &
      gzip ../${IDENTIFIER}_cat_R2.fastq &
      wait

      echo "Attempting to DELETE old files"
      if [[ $IDENTIFIER == "" ]]; then
          echo "ISSUE deleting files (identifier not set)"
          exit
      else
          if [[ (-f "../${IDENTIFIER}_cat_R1.fastq.gz") && (-f "../${IDENTIFIER}_cat_R2.fastq.gz") ]]; then
            echo "Prodeeding with DELETING old fastq files"
            rm $datapath/${IDENTIFIER}/*.fastq
          else
            echo "ISSUE deleting files (.gz files not found)"
            exit
          fi
      fi
        
      echo "${IDENTIFIER} job done."
      echo "${IDENTIFIER} job done." > ../log2_${IDENTIFIER}.txt
        
done










