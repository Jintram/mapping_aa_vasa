#!/bin/bash

datapath='/hpc/hub_oudenaarden/mwehrens/fastq/WANG4/'
#datapath='/hpc/hub_oudenaarden/mwehrens/fastq/WANG/'
#datapath=/Volumes/workdrive_m.wehrens_hubrecht/Fastq__raw_files/Wang/

for IDENTIFIER in p.N1.plate.97474.part.1 p.N1.plate.97474.part.2 p.N2.plate.97452.part.1 p.N2.plate.97452.part.2 p.N2.plate.97493.part.1 p.N2.plate.97493.part.2 p.N3.plate.97438.part.1 p.N3.plate.97438.part.2 p.N4.plate.97461.part.1 p.N4.plate.97461.part.2 p.N5.plate.97458.part.1 p.N5.plate.97458.part.2 p.N13.plate.100355.part.1 p.N13.plate.100355.part.2 p.N14.plate.104720.part.1 p.N14.plate.104720.part.2
#for IDENTIFIER in GSM2970361_N1_LV GSM2970358_N2_LV GSM2970362_N3_LV GSM2970366_N4_LV GSM2970360_N5_LV GSM2970359_N6_LA GSM2970369_N7_LA GSM2970368_N8_LA GSM2970365_N9_LA GSM2970363_N10_LA GSM2970364_N11_LA GSM2970367_N12_LA GSM3449619_N13 GSM3449620_N14
do

      check_FLAG=1

      cd $datapath/${IDENTIFIER}/

      filename=${IDENTIFIER}_SRR_Acc_List.txt

      echo "Checking job ${IDENTIFIER} .."
      echo "Checking job ${IDENTIFIER} .." > check_${IDENTIFIER}.txt

      n=1
      while read line; do
        
          # check read 1 file
          subfile1_to_check="${line}_1.fastq"
          if [ -f "${subfile1_to_check}" ]; then
              echo "$subfile1_to_check exists." >> check_${IDENTIFIER}.txt
          else 
              echo "$subfile1_to_check is missing!"
              echo "$subfile1_to_check does not exist." >> check_${IDENTIFIER}.txt
              check_FLAG=0
          fi
          
          # check read 2 file
          subfile2_to_check="${line}_2.fastq"
          if [ -f "${subfile2_to_check}" ]; then
              echo "$subfile2_to_check exists." >> check_${IDENTIFIER}.txt
          else 
              echo "$subfile2_to_check is missing!"
              echo "$subfile2_to_check does not exist." >> check_${IDENTIFIER}.txt
              check_FLAG=0
          fi
        
      done < $filename

      if [[ ${check_FLAG} == 1 ]]; then
          echo "Check PASSED."
          echo "Check PASSED." >> check_${IDENTIFIER}.txt
      else
          echo "Check NOT passed."
          echo "Check NOT passed." >> check_${IDENTIFIER}.txt 
      fi

done










