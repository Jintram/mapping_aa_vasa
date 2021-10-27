#!/bin/bash

if [ $# -ne 4 ]
then
    #echo "To L_TS_deal_with_singleandmultimappers_paired_stranded.sh, please, give:"
    #echo "1) input bam file"
    #echo "2) bed file for introns, exons and tRNA"
    #echo "3) stranded protocol (y/n)"
    #echo "4) paired data (y/n)"
    echo "Please, give (1) general param file, (2) run param file, (3) input bam file, (4) paired [y/n]"    
    exit
fi



#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

################################################################################

general_parameter_filepath=$1
run_parameter_filepath=$2
inbam=$3
paired=$4

source $general_parameter_filepath
source $run_parameter_filepath
current_dir=$(pwd)

if [[ $notempdir == "" ]]; then # corrected
  echo "Going to $TMPDIR"
  cd $TMPDIR
else
  echo "Going to $outdir"
  cd $outdir
fi

################################################################################
# Outline of this script
# 1) Extract cigar and nM (mismatches) parameter (add it to name entry)
# 2) Convert to bed
# 3) Calculate feature overlap between genome annotation file
# 4) Annotate extra properties of those overlaps (like does it fall completely in a feature, etc.)
#
# This version doesn't have "non-stranded" support (Anna's version did)
################################################################################

################################################################################
# (1) Extract cigar and nM (mismatches) parameter (add it to name entry)
# Updated version of adding extra info by MW in python

if [ $paired == 'y' ]
then
  ${p2samtools}/samtools view -h -f 2 -o ${inbam%.bam}.f2.sam $inbam
  f2_str='.f2'
  # Important!: the "-f 2" option removes singletons from the mapping;
  # these can occur because STAR will map one mate only if other mate's length
  # is <1/3 that one.
  # See also: https://groups.google.com/g/rna-star/c/K8yVdkTlWoY
else
  # if not paired, don't use "-f 2" option 
  ${p2samtools}/samtools view -h -o ${inbam%.bam}.sam $inbam
  f2_str=''
fi

# The following python code can split R1 and R2 into separate files, or 
# just throw everything in one file
# Split in two files:
# $pythonbin ${p2s}/TS_mw_extract_cigar_nM.py ${inbam%.bam}.f2.sam ${inbam%.bam}.f2.singlemappers.sam ${inbam%.bam}.f2.multimappers.sam '.R1;.R2'
# Everything in one file
$pythonbin ${p2s}/TS_mw_extract_cigar_nM.py ${inbam%.bam}${f2_str}.sam ${inbam%.bam}${f2_str}.singlemappers.sam ${inbam%.bam}${f2_str}.multimappers.sam

# Convert back to bam
${p2samtools}/samtools view -S -b -o ${inbam%.bam}${f2_str}.singlemappers.bam ${inbam%.bam}${f2_str}.singlemappers.sam 
${p2samtools}/samtools view -S -b -o ${inbam%.bam}${f2_str}.multimappers.bam ${inbam%.bam}${f2_str}.multimappers.sam

################################################################################
# Previous version by Anna (didn't deal with paired end information)
#samtools view -h $inbam | awk 'BEGIN{OFS="\t"} {
#    if ($1 ~ /^@/) {print $0}                    # ignore but print first part of bam file
#  else if ($0 ~ /NH:i:1\tHI:i:1\t/) {          # if single mapper (defined by NH:i:1)
#        for (i=1; i<=NF; i++) {                  # loop over entries
#            if ($i ~ /nM:i:[0-9]/) {             # collect nM stats (# mismatches)
#                col=i; nm=substr($col, 6, length($col)) # ignore "nM:i:" part and take only the number
#            }
#        };
#        $1=$1";CG:"$6";nM:"nm; print $0          # now output line again, but add info to read name
#    }
#}' | samtools view -Sb > ${inbam%.bam}.singlemappers.bam

# TODO: CHECK WHETHER IT'S NECESSARY TO EDIT THIS CODE TO DEAL WITH PAIRED END READS
#if [ $paired == "y" ]
#then
#  samtools sort -n -o ${inbam%.bam}.singlemappers.nsorted.bam ${inbam%.bam}.singlemappers.bam
#  $p2b/pairToBed -abam ${inbam%.bam}.singlemappers.nsorted.bam -b $refBED > test.txt # | $p2b/bedtools sort > ${inbam%.bam}.singlemappers.bed
#else

################################################################################
# (2) Convert to bed (this can also handle pairs, which will be put together per set of two mates)

# Strategy I: sorted by index (Anna did this, necessary?)
#$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.singlemappers.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.singlemappers.bed
#$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.multimappers.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.multimappers.bed

# Strategy II: sorted by pairs:
$p2b/bedtools bamtobed -i ${inbam%.bam}${f2_str}.singlemappers.bam > ${inbam%.bam}${f2_str}.singlemappers.bed
$p2b/bedtools bamtobed -i ${inbam%.bam}${f2_str}.multimappers.bam > ${inbam%.bam}${f2_str}.multimappers.bed

# Strategy III: What exactly does pairtobed do?
# $p2b/pairToBed -bedpe -abam ${inbam%.bam}.f2.singlemappers.bam -b $refBED > ${inbam%.bam}.f2.singlemappers.pair.bed # | $p2b/bedtools sort > ${inbam%.bam}.singlemappers.bed

# Strategy IV: do it with reads R1 and R2 being separated into separate files (see above)
# Maybe it's easier to just ignore all of this crap, and do the files manually, separately
#samtools view -S -b -o ${inbam%.bam}.f2.singlemappers.R1.bam ${inbam%.bam}.f2.singlemappers.R1.sam 
#samtools view -S -b -o ${inbam%.bam}.f2.singlemappers.R2.bam ${inbam%.bam}.f2.singlemappers.R2.sam 
#samtools view -S -b -o ${inbam%.bam}.f2.multimappers.R1.bam ${inbam%.bam}.f2.multimappers.R1.sam
#samtools view -S -b -o ${inbam%.bam}.f2.multimappers.R2.bam ${inbam%.bam}.f2.multimappers.R2.sam

# Strategy V: separate files, but including sorting
# First for single mappers
#$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.singlemappers.R1.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.singlemappers.R1.bed
#$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.singlemappers.R2.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.singlemappers.R2.bed
# Now for the multi-mappers
#$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.multimappers.R1.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.multimappers.R1.bed
#$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.multimappers.R2.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.multimappers.R2.bed

################################################################################
# (3) Then calculate intersect

# Strategy IV: (mates split)
#$p2b/bedtools intersect -a ${inbam%.bam}.f2.singlemappers.R1.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.singlemappers.R1.intersect.bed
#$p2b/bedtools intersect -a ${inbam%.bam}.f2.singlemappers.R2.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.singlemappers.R2.intersect.bed
#$p2b/bedtools intersect -a ${inbam%.bam}.f2.multimappers.R1.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.multimappers.R1.intersect.bed
#$p2b/bedtools intersect -a ${inbam%.bam}.f2.multimappers.R2.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.multimappers.R2.intersect.bed


# what happens if we do this for the paired bed file? # See here XXX!!!!!!!
$p2b/bedtools intersect -a ${inbam%.bam}${f2_str}.singlemappers.bed -b $refBED -wa -wb > ${inbam%.bam}${f2_str}.singlemappers.intersect.bed
$p2b/bedtools intersect -a ${inbam%.bam}${f2_str}.multimappers.bed -b $refBED -wa -wb > ${inbam%.bam}${f2_str}.multimappers.intersect.bed
  # output will be (multiple) features, sorted by feature overlap of the reads
  # mates are labeled by /1 and /2 at the end of the read name
  # note that IGV doesn't have the comprehensive set of genes,
  # so to retrieve some features, look at http://www.ensembl.org/
  
################################################################################

# Now update each feature to add information about how it overlaps with features
# Note that the x-parameter is strand-specific
# updated with paired-cigar appropriate search term "/;CG:|;C1:/"
# updated with R1/R2 detection and appropriate handling of strands
if [[ $stranded == 'y' ]]; then
  strandedyesno_code_part1='if ((readstrand==refstrand && XM==2)||(readstrand!=refstrand && XM==1)) {'
  strandedyesno_code_part2='}'
else
  strandedyesno_code_part1=''
  strandedyesno_code_part2=''
fi

awk_pt1='BEGIN {OFS="\t"; w="T"} {
    # annotate reads to know fall inbetween inton/exon boundaries, etc
    
    # jS = junction
    # 5 & 3 --> which end 
    # IN:intron
    # Cigar has splicing information
    # w: write to file (True/False)
  
    # give fields names
    chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;

    # obtain whether R1 or R2
    # (mw) first obtain the /1 or /2 part (if there), which becomes part of the bed files for paired reads (also remove that part)
    if (readname ~ /\/1$|\/2$/) { 
      sidx=match(readname, /\/1$|\/2$/); XM=substr(readname,sidx+1,1); readname=substr(readname,0,sidx-1) 
    } else { XM = 2 } # when theres no annotation for R1/R2, it is assumed were dealing with Read 2.
    
    # separate readname into two parts (leaving nM/NM into part 1, CG and other fields into part 2)
    sx=match(readname, /;CG:|;C1:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
    '
   
awk_pt2=$strandedyesno_code_part1       
   
awk_pt3='
        # then check whether read fall entirely into annotation
        # if so, add "IN" to signal that
        if ((readstart >= refstart) && (readend <= refend)) {
            readname=readname";jS:IN"; w="T"
        # conversely, if both ends fall outside ref, annotate that, but in fact dont even write it (ie filter it out)
        } else if ((readstart < refstart) && (readend > refend)) {
                readname=readname";jS:OUT"; w="F";
        # if one side falls outside, annotate which part falls out
        } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
                    if (readstrand=="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
                } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
                    if (readstrand=="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
                } else {print $0 > "checkme.txt"
        }
    
        # calculate x parameter (relative position?)
        if (readstrand=="+") {
          x=1-(geneend-readend)/genelen
        } else {
          x=1-(readstart-genestart)/genelen
        }
        
        # again separate the readname into two parts
        sx=match(readname, /;CG:|;C1:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))";XM:"XM
        
        if (w=="T") {
          print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
        }
      '
awk_pt4=$strandedyesno_code_part2
    
final_awk_command=$(echo "${awk_pt1}${awk_pt2}${awk_pt3}${awk_pt4}")

echo "************************************************************"
echo "(Debugging, remove..)"
echo "FINAL AWK COMMAND:"
echo "************************************************************"
echo "${final_awk_command}"
echo "************************************************************"
    
    
#awk_command_singlemappers_part2_R2='    
#    # only look at feature when annotation is on same strand as mapped read
#    if (readstrand==refstrand) {'
#awk_command_singlemappers_part2_R1='    
#        # only look at feature when annotation is on opposite strand as mapped read (R1 is on opposite of original gene strand)
#        if (readstrand!=refstrand) {'       
    

################################################################################

# Note that the bedtools "intersect" option is used to obtain overlap between
# the annotation reference and mapped reads.
# This is processed further using awk, and then other scripts.
  #cat ${inbam%.bam}.f2.singlemappers.intersect.bed | awk $awk_command_singlemappers > ${inbam%.bam}.singlemappers_genes.bed

cat ${inbam%.bam}${f2_str}.singlemappers.intersect.bed | awk "${final_awk_command}" > ${inbam%.bam}.singlemappers_genes.bed
#sort -k4 ${inbam%.bam}.singlemappers_genes.bed > ${inbam%.bam}.nsorted.singlemappers_genes.bed 
  # sorting here is not necessary

cat ${inbam%.bam}${f2_str}.multimappers.intersect.bed | awk "${final_awk_command}" > ${inbam%.bam}.multimappers_genes.bed
sort -k4 ${inbam%.bam}.multimappers_genes.bed > ${inbam%.bam}.nsorted.multimappers_genes.bed 
  # note that the python script to process these reads is completely dependent on reads
  # being sorted by name, if not, reads will be double counted

################################################################################

if [[ -f ${inbam%.bam}.singlemappers_genes.bed && -f ${inbam%.bam}.nsorted.multimappers_genes.bed ]]; then
  if [[ $nocleanup = "" ]]; then
    echo "cleaning up some files" # prevent this by setting "nocleanup"
    
    #rm ${inbam%.bam}.singlemappers_genes.bed
    rm ${inbam%.bam}.multimappers_genes.bed
    
    rm $inbam
    rm ${inbam%.bam}${f2_str}.singlemappers.bam
    rm ${inbam%.bam}${f2_str}.multimappers.bam
    
    rm ${inbam%.bam}${f2_str}.sam
    rm ${inbam%.bam}${f2_str}.singlemappers.sam
    rm ${inbam%.bam}${f2_str}.multimappers.sam
  fi
else
  echo "output files not detected; returning non-zero exit"
  exit 1  
fi

################################################################################






# OLD VERSION
# First singlemappers  
#cat ${inbam%.bam}.f2.singlemappers.R1.intersect.bed | awk $awk_command_singlemappers_R1 > ${inbam%.bam}.R1.singlemappers_genes.bed
#cat ${inbam%.bam}.f2.singlemappers.R2.intersect.bed | awk $awk_command_singlemappers_R2 > ${inbam%.bam}.R2.singlemappers_genes.bed
#sort -k4 ${inbam%.bam}.R1.singlemappers_genes.bed > ${inbam%.bam}.R1.nsorted.singlemappers_genes.bed
#sort -k4 ${inbam%.bam}.R2.singlemappers_genes.bed > ${inbam%.bam}.R2.nsorted.singlemappers_genes.bed
# Then multimappers
#cat ${inbam%.bam}.f2.multimappers.R1.intersect.bed | awk $awk_command_singlemappers_R1 > ${inbam%.bam}.R1.multimappers_genes.bed
#cat ${inbam%.bam}.f2.multimappers.R2.intersect.bed | awk $awk_command_singlemappers_R2 > ${inbam%.bam}.R2.multimappers_genes.bed
#sort -k4 ${inbam%.bam}.R1.multimappers_genes.bed > ${inbam%.bam}.R1.nsorted.multimappers_genes.bed
#sort -k4 ${inbam%.bam}.R2.multimappers_genes.bed > ${inbam%.bam}.R2.nsorted.multimappers_genes.bed

# some tests (note that numbers will also depend on # features within a gene that happen to overlap)
#grep "_MYBPC3_" ${inbam%.bam}.R1.multimappers_genes.bed | wc -l
# 0
#grep "_MYBPC3_" ${inbam%.bam}.R2.multimappers_genes.bed | wc -l
# 6
#grep "_MYBPC3_" ${inbam%.bam}.R1.singlemappers_genes.bed | wc -l
# 127
#grep "_MYBPC3_" ${inbam%.bam}.R2.singlemappers_genes.bed | wc -l
# 188
# So why the discrepancy? For singlemappers, probably because 
# R2 is longer read and thus has more features..
# But for the multimappers, it's weird that an MYBPC3 R2 is not matched
# by an R1 MYBPC3.

# For GAPDH, it seems reasonably consistent
#grep "_GAPDH_" ${inbam%.bam}.R1.multimappers_genes.bed | wc -l
# 35
#grep "_GAPDH_" ${inbam%.bam}.R2.multimappers_genes.bed | wc -l
# 37
#grep "_GAPDH_" ${inbam%.bam}.R1.singlemappers_genes.bed | wc -l
# 2014
#grep "_GAPDH_" ${inbam%.bam}.R2.singlemappers_genes.bed | wc -l
# 2032


#sort -k4 ${inbam%.bam}.singlemappers_genes.bed > ${inbam%.bam}.nsorted.singlemappers_genes.bed # unnecessary

