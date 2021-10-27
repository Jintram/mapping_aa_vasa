#!/bin/sh

# Install necessary scripts using:
#
# conda install -c bioconda subread
# conda install -c bioconda umi_tools
# @HPC, umi_tools needed to be installed with pip
# pip install umi-tools

# This follows
# https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html

example="
#bamfile=/Volumes/fastq_m.wehrens/Mapping/WANG4/mapping/head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out
#gtffile=/Volumes/fastq_m.wehrens/Mapping/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
bamfile=/hpc/hub_oudenaarden/mwehrens/fastq/WANG4/mapping/head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out
nrthreads=1
trimUMI=0
steptorun=1
countmulti=no
"

# To run, e.g. execute:
# myscript=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/countTables_umiTools.sh
# myscript=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04/countTables_umiTools.sh
# $myscript $bamfile $gtffile $nrthreads $trimUMI 1
# $myscript $bamfile $gtffile $nrthreads 10 2

if [[ $bamfile == "" ]]; then
  bamfile=$1
  gtffile=$2
  nrthreads=$3
  trimUMI=$4
  steptorun=$5
  countmulti=$6
fi

echo "Going ($steptorun).."

# Testing the UMI tools thing
# conda install -c bioconda subread

if [[ $steptorun == "1" ]]; then

    echo "Running step 1"

    # Run featureCounts
    if [[ $countmulti == "yes" ]]; then
      # see "featureCounts --help" for info, multis are counted towards all hits
      multioption="-M"
    elif [[ $countmulti == "no" ]]; then
      # default behavior is to ignore multi mappers, 
      # so no additional option needed
      multioption="" 
    else
      exit 1
    fi
      
    featureCounts -a ${gtffile} \
                  -o gene_assigned \
                  -R BAM ${bamfile}.bam \
                  -T ${nrthreads} \
                  -s 1 ${multioption}
        
fi                  
        
if [[ $steptorun == "2" ]]; then
          
    echo "Running step 2"        
            
    ################################################################################
    # Addition by MW

    awk_command_rearrange='BEGIN{FS="\t"} {
      if (substr($0,1,1)=="@") {print $0} 
      else {
        split($1,first_part,";"); 
        printf first_part[1] "\t"; 
        for (i=2; i<=NF; i++) 
          {printf $i "\t"}; 
        for (i=2; i<length(first_part); i++) 
          {gsub(/:/, ":Z:", first_part[i]); printf first_part[i] "\t"}; 
        gsub(/:/, ":Z:", first_part[length(first_part)]); print first_part[length(first_part)]
      }}'
    awk_command_trimUMI='BEGIN{FS="\t"; OFS="\t"; UMI_TRIM_LN=10} {
        for (i=1; i<=NF; i++) {
            if(substr($i,1,2)=="RX") {$i=substr($i,1,UMI_TRIM_LN+5)}        
        };
        print $0
      }'  

    ################################################################################
    # Unpiped version

    #      # TODO: perhaps this should be one big pipe!!
    #      # First create a sam file
    #      samtools view -h -o ${bamfile}.bam.featureCounts.sam ${bamfile}.bam.featureCounts.bam    
    #      # We need to extract the tags that we added to the title and add them as separate field
    #      # awk 'BEGIN{FS="\t"} {if (substr($0,1,1)=="@") {print $0} else {split($1,first_part,";"); printf first_part[1] "\t"; for (i=2; i<=NF; i++) {printf $i "\t"}; for (i=2; i<length(first_part); i++) {printf first_part[i] "\t"}; print first_part[length(first_part)]}}' ${bamfile}.bam.featureCounts.sam > ${bamfile}.bam.featureCounts_MW.sam
    #      # Re-arrange the bam file
    #      awk $awk_command_rearrange ${bamfile}.bam.featureCounts.sam > ${bamfile}.bam.featureCounts_MW.sam
    #      # Perform UMI trimming if necessary
    #      awk $awk_command_trimUMI ${bamfile}.bam.featureCounts_MW.sam > ${bamfile}.bam.featureCounts_MW_trimmed.sam
    #      mv ${bamfile}.bam.featureCounts_MW_trimmed.sam ${bamfile}.bam.featureCounts_MW.sam
    #      # Convert this back to bam
    #      samtools view -S -b ${bamfile}.bam.featureCounts_MW.sam > ${bamfile}.bam.featureCounts_MW.bam

    ################################################################################
    # Piped version of above

    # Rearrange the bam file to accomodate UMI tools input
    samtools view -h ${bamfile}.bam.featureCounts.bam | awk "$awk_command_rearrange" | samtools view -S -b > ${bamfile}.bam.featureCounts_MW.bam

    # Trim UMIs if desired (this is awkward solution to UMIs not being equal length)
    # one could also deliberate running umi-tools without the python assert option..
    if [[ $trimUMI > 0 ]]; then
      echo "Trimming UMIs to deal with unequal length UMIs.."
      samtools view -h ${bamfile}.bam.featureCounts_MW.bam | awk "$awk_command_trimUMI" | samtools view -S -b > ${bamfile}.bam.featureCounts_MW_umitrimmed.bam
      mv ${bamfile}.bam.featureCounts_MW_umitrimmed.bam ${bamfile}.bam.featureCounts_MW.bam  
    fi

    ################################################################################
                    
    # Sorting required for UMI tools
    samtools sort ${bamfile}.bam.featureCounts_MW.bam -o ${bamfile}.assigned_sorted_MW.bam
    samtools index ${bamfile}.assigned_sorted_MW.bam

    # Original command from tutorial
    # umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${bamfile}.assigned_sorted.bam -S ${bamfile}.counts.tsv.gz

    # Now produce the count tables.
    umi_tools count --cell-tag=SM --extract-umi-method=tag --umi-tag=RX --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${bamfile}.assigned_sorted_MW.bam -S ${bamfile}.counts.tsv.gz

fi

if [[ $steptorun == "3" ]]; then
  
    echo "Step 3: cleaning up some temp files."
  
    # only if final output exists 
    if [ -f "${bamfile}.counts.tsv.gz" ]; then
      echo "Removing temporary file ${bamfile}.bam.featureCounts_MW.bam."
      rm ${bamfile}.bam.featureCounts_MW.bam
    fi

fi

echo "End of script reached."





