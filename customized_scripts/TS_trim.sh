#!/bin/bash

if [ $# -ne 3 ]
then
    #echo "Please, give (1) input fastq file; (2) path2trimgalore; (3) path2cutadapt;"    
    echo "Please, give (1) general param file, (2) run param file, (3) file to trim"    
  exit
fi

#file2trim=$1
#path2trimgalore=$2
#path2cutadapt=$3

general_parameter_filepath=$1
run_parameter_filepath=$2
file2trim=$3

source $general_parameter_filepath
source $run_parameter_filepath
current_dir=$(pwd)
echo "Going to $outdir"
cd $outdir

# trim adaptors
${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt ${file2trim}_cbc.fastq.gz
mv ${file2trim}_cbc.fastq.gz_trimming_report.txt ${file2trim}_trimming_report.txt

# trim homopolymers 
${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a AAAAAAAAAAAAAAAAAAAA --no_report_file ${file2trim}_cbc_trimmed.fq.gz
mv ${file2trim}_cbc_trimmed_trimmed.fq.gz ${file2trim}_cbc_trimmed_HA.fq.gz # renaming output file
rm ${file2trim}_cbc_trimmed.fq.gz # removal of input file

${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a TTTTTTTTTTTTTTTTTTTT --no_report_file ${file2trim}_cbc_trimmed_HA.fq.gz
mv ${file2trim}_cbc_trimmed_HA_trimmed.fq.gz ${file2trim}_cbc_trimmed_HAT.fq.gz
rm ${file2trim}_cbc_trimmed_HA.fq.gz

${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a CCCCCCCCCCCCCCCCCCCC --no_report_file ${file2trim}_cbc_trimmed_HAT.fq.gz
mv ${file2trim}_cbc_trimmed_HAT_trimmed.fq.gz ${file2trim}_cbc_trimmed_HATC.fq.gz
rm ${file2trim}_cbc_trimmed_HAT.fq.gz

${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt -a GGGGGGGGGGGGGGGGGGGG --no_report_file ${file2trim}_cbc_trimmed_HATC.fq.gz
mv ${file2trim}_cbc_trimmed_HATC_trimmed.fq.gz ${file2trim}_cbc_trimmed_HATCG.fq.gz
rm ${file2trim}_cbc_trimmed_HATC.fq.gz

# File cleanup
if [[ $nocleanup = "" ]]; then
  echo "cleaning up some files" # prevent this by setting "nocleanup"
  rm ${file2trim}_cbc.fastq.gz
fi

cd $current_dir



#java -jar ${path2trimmomatic}/trimmomatic-0.36.jar SE -phred33 ${file2trim%.fastq.gz}_trimmed.fq.gz ${file2trim%.fastq.gz}_trimHomo.fq ILLUMINACLIP:/hpc/hub_oudenaarden/fsalmen/barcodes/homopolymers.fa:0::0:0 MINLEN:20

#gzip ${file2trim%.fastq.gz}_trimHomo.fq


