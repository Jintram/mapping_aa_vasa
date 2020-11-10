#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) input bam file"
    echo "2) bed file for introns, exons and tRNA"
    echo "3) stranded protocol (n/y)"
    echo "4) paired (n/y)"
    exit
fi

inbam=$1
refBED=$2
stranded=$3
paired=$4
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

#############
# Updated version of adding extra info by MW in python

samtools view -h -f 2 -o ${inbam%.bam}.f2.sam $inbam
  # Important!: the "-f 2" option removes singletons from the mapping;
  # these can occur because STAR will map one mate only if other mate's length
  # is <1/3 that one.
  # See also: https://groups.google.com/g/rna-star/c/K8yVdkTlWoY
$pythonbin ${p2s}/TS_mw_extract_cigar_nM.py ${inbam%.bam}.f2.sam ${inbam%.bam}.f2.singlemappers.sam ${inbam%.bam}.f2.multimappers.sam '.R1;.R2'

# Convert back to bam
samtools view -S -b -o ${inbam%.bam}.f2.singlemappers.bam ${inbam%.bam}.f2.singlemappers.sam 
samtools view -S -b -o ${inbam%.bam}.f2.multimappers.bam ${inbam%.bam}.f2.multimappers.sam

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

# For local use, we also need to give path
# On server, probably "bedtools ..." will work as command
p2b=/Users/m.wehrens/Software_custom/bedtools2/bin/

# TODO: CHECK WHETHER IT'S NECESSARY TO EDIT THIS CODE TO DEAL WITH PAIRED END READS
#if [ $paired == "y" ]
#then
#  samtools sort -n -o ${inbam%.bam}.singlemappers.nsorted.bam ${inbam%.bam}.singlemappers.bam
#  $p2b/pairToBed -abam ${inbam%.bam}.singlemappers.nsorted.bam -b $refBED > test.txt # | $p2b/bedtools sort > ${inbam%.bam}.singlemappers.bed
#else
$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.singlemappers.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.singlemappers.bed
#fi

######################
# TESTING OUT SOME STUFF

# Maybe we don't want to sort this actually, because the reads are grouped per read pair
$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.singlemappers.bam > ${inbam%.bam}.f2.singlemappers.bed

# What exactly does pairtobed do?
$p2b/pairToBed -bedpe -abam ${inbam%.bam}.f2.singlemappers.bam -b $refBED > ${inbam%.bam}.f2.singlemappers.pair.bed # | $p2b/bedtools sort > ${inbam%.bam}.singlemappers.bed
######################

# Maybe it's easier to just ignore all of this crap, and do the files manually, separately
samtools view -S -b -o ${inbam%.bam}.f2.singlemappers.R1.bam ${inbam%.bam}.f2.singlemappers.R1.sam 
samtools view -S -b -o ${inbam%.bam}.f2.singlemappers.R2.bam ${inbam%.bam}.f2.singlemappers.R2.sam 
samtools view -S -b -o ${inbam%.bam}.f2.multimappers.R1.bam ${inbam%.bam}.f2.multimappers.R1.sam
samtools view -S -b -o ${inbam%.bam}.f2.multimappers.R2.bam ${inbam%.bam}.f2.multimappers.R2.sam

# First for single mappers
# Convert to bed and sort
$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.singlemappers.R1.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.singlemappers.R1.bed
$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.singlemappers.R2.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.singlemappers.R2.bed
# Then calculate intersect
$p2b/bedtools intersect -a ${inbam%.bam}.f2.singlemappers.R1.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.singlemappers.R1.intersect.bed
$p2b/bedtools intersect -a ${inbam%.bam}.f2.singlemappers.R2.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.singlemappers.R2.intersect.bed

# Now for the multi-mappers
# Convert to bed and sort
$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.multimappers.R1.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.multimappers.R1.bed
$p2b/bedtools bamtobed -i ${inbam%.bam}.f2.multimappers.R2.bam | $p2b/bedtools sort > ${inbam%.bam}.f2.multimappers.R2.bed
# Then calculate intersect
$p2b/bedtools intersect -a ${inbam%.bam}.f2.multimappers.R1.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.multimappers.R1.intersect.bed
$p2b/bedtools intersect -a ${inbam%.bam}.f2.multimappers.R2.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.multimappers.R2.intersect.bed


# what happens if we do this for the paired bed file? # See here XXX!!!!!!!
See here XXX!!!!!!!
$p2b/bedtools intersect -a ${inbam%.bam}.f2.singlemappers.bed -b $refBED -wa -wb > ${inbam%.bam}.f2.singlemappers.intersect.bed
  # i don't quite get how it sorts the reads...
  
################################################################################

# updated with paired-cigar appropriate search term "/;CG:|;C1:/"
awk_command_singlemappers_part1='BEGIN {OFS="\t"; w="T"} {
    # annotate reads to know fall inbetween inton/exon boundaries, etc
    
    # jS = junction
    # 5 & 3 --> which end 
    # IN:intron
    # Cigar has splicing information
    # w: write to file (True/False)
  
    # give fields names
    chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;
    
    # (mw) first remove the /1 or /2 part (if there), which becomes part of the bed files for paired reads
    # it isnt applicable since we split R1 and R2 in two separate files
    if (readname ~ /\/1$|\/2$/) { sidx=match(readname, /\/1$|\/2$/); readname=substr(readname,0,sidx-1) }
    
    # separate readname into two parts (leaving nM/NM into part 1, CG and other fields into part 2)
    sx=match(readname, /;CG:|;C1:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))'
awk_command_singlemappers_part2_R2='    
    # only look at feature when annotation is on same strand as mapped read
    if (readstrand==refstrand) {'
awk_command_singlemappers_part2_R1='    
        # only look at feature when annotation is on opposite strand as mapped read (R1 is on opposite of original gene strand)
        if (readstrand!=refstrand) {'    
awk_command_singlemappers_part3='      
        # then check whether read fall entirely into annotation
        # if so, add "IN" to signal that
        if ((readstart >= refstart) && (readend <= refend)) {
            readname=readname";jS:IN"; w="T"
        # conversely, if both ends fall outside ref, annotate that, but in fact dont even write it (ie filter it out)
        } else if ((readstart < refstart) && (readend > refend)) {
                readname=readname";jS:OUT"; w="F";
        # if one side falls outside, annotate which part falls out
        } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
                    if (readstrand="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
                } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
                    if (readstrand="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
                } else {print $0 > "checkme.txt"
        }
    
        # calculate x parameter
        if (readstrand=="+") {
          x=1-(geneend-readend)/genelen
        } else {
          x=1-(readstart-genestart)/genelen
        }
        
        # again separate the readname into two parts
        sx=match(readname, /;CG:|;C1:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        
        if (w=="T") {
          print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
        }
    }'
# difference between the two lies in the strands
awk_command_singlemappers_R1=${awk_command_singlemappers_part1}${awk_command_singlemappers_part2_R1}${awk_command_singlemappers_part3}
awk_command_singlemappers_R2=${awk_command_singlemappers_part1}${awk_command_singlemappers_part2_R2}${awk_command_singlemappers_part3}

################################################################################

# Note that the bedtools "intersect" option is used to obtain overlap between
# the annotation reference and mapped reads.
# This is processed further using awk, and then other scripts.
  #cat ${inbam%.bam}.f2.singlemappers.intersect.bed | awk $awk_command_singlemappers > ${inbam%.bam}.singlemappers_genes.bed

# First singlemappers  
cat ${inbam%.bam}.f2.singlemappers.R1.intersect.bed | awk $awk_command_singlemappers_R1 > ${inbam%.bam}.R1.singlemappers_genes.bed
cat ${inbam%.bam}.f2.singlemappers.R2.intersect.bed | awk $awk_command_singlemappers_R2 > ${inbam%.bam}.R2.singlemappers_genes.bed
sort -k4 ${inbam%.bam}.R1.singlemappers_genes.bed > ${inbam%.bam}.R1.nsorted.singlemappers_genes.bed
sort -k4 ${inbam%.bam}.R2.singlemappers_genes.bed > ${inbam%.bam}.R2.nsorted.singlemappers_genes.bed

# Then multimappers
cat ${inbam%.bam}.f2.multimappers.R1.intersect.bed | awk $awk_command_singlemappers_R1 > ${inbam%.bam}.R1.multimappers_genes.bed
cat ${inbam%.bam}.f2.multimappers.R2.intersect.bed | awk $awk_command_singlemappers_R2 > ${inbam%.bam}.R2.multimappers_genes.bed
sort -k4 ${inbam%.bam}.R1.multimappers_genes.bed > ${inbam%.bam}.R1.nsorted.multimappers_genes.bed
sort -k4 ${inbam%.bam}.R2.multimappers_genes.bed > ${inbam%.bam}.R2.nsorted.multimappers_genes.bed

# some tests 
grep "_MYBPC3_" ${inbam%.bam}.R1.multimappers_genes.bed | wc -l
# 0
grep "_MYBPC3_" ${inbam%.bam}.R2.multimappers_genes.bed | wc -l
# 6
grep "_MYBPC3_" ${inbam%.bam}.R1.singlemappers_genes.bed | wc -l
# 127
grep "_MYBPC3_" ${inbam%.bam}.R2.singlemappers_genes.bed | wc -l
# 188
# So why the discrepancy? For singlemappers, probably because 
# R2 is longer read and thus has more features..
# But for the multimappers, it's weird that an MYBPC3 R2 is not matched
# by an R1 MYBPC3.

# For GAPDH, it seems reasonably consistent
grep "_GAPDH_" ${inbam%.bam}.R1.multimappers_genes.bed | wc -l
# 35
grep "_GAPDH_" ${inbam%.bam}.R2.multimappers_genes.bed | wc -l
# 37
grep "_GAPDH_" ${inbam%.bam}.R1.singlemappers_genes.bed | wc -l
# 2014
grep "_GAPDH_" ${inbam%.bam}.R2.singlemappers_genes.bed | wc -l
# 2032


#sort -k4 ${inbam%.bam}.singlemappers_genes.bed > ${inbam%.bam}.nsorted.singlemappers_genes.bed # unnecessary

