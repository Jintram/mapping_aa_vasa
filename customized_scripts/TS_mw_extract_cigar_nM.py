#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:26:36 2020

@author: m.wehrens
"""

# This files extracts the cigar string and number of misreads,
# and adds it to the read name, to be used later.
#
# Currently it also does quite some extra checks, since I was
# debugging; TODO: remove those.

# Notes:
# STAR outputs nM and NM, from the manual (v 2.7.3a, p.11):
# nM : is the number of mismatches per (paired) alignment, not to be confused with NM, which is
# the number of mismatches in each mate.
# 
# I read somewhere that STAR might also output a combined CIGAR string, 
# but that doesn't seem to be the case here. So I'll just combine 
# those with a slash.. (Lp was suggested format, where L is a number representing the gap between pairs -
# https://www.biostars.org/p/113136/)
#
# NH = Number of reported alignments that contain the query in the current record
#
# From STAR manual:
# --outSAMtype BAM Unsorted
# output unsorted Aligned.out.bam The paired ends of an alignment are always adjacent,
# and multiple alignments of a read are adjacent as well. This "unsorted" can be directly
# used with downstream software such as HTseq, without the need of name sorting. The order
# of the reads will match that of the input FASTQ(A)s only if one thread is used
#
# Note that STAR will output reads where only one of the mates is mapped as singletons
# if one read size < 1/3 of the other.
# These singletons can be filtered out beforehand using the "-f 2" option in samtools.
#
# ===
# Excerpt from Google groups:
# STAR controls the mapped length/score of the output reads, requiring >66% of the total read length (sum of mates) mapped.
# This prevents output of unpaired alignments if the lengths of the mates are similar. However, is the length of one mate is <1/3 of the other, this filter will allow singletons.
# To make sure that no singletons are output, please make --outFilterMatchNmin <MaxMateLength+1>
#
# If you do not mind running samtoold, a better solution is to filter the singletons after mapping:
# $ samtools view -b -f 0x2 Aligned.out.bam > Aligned.PEonly.bam
# This will work for both unsorted and sorted BAMs.
# 
#                   -- https://groups.google.com/g/rna-star/c/K8yVdkTlWoY
# ===

import re, sys

try:
    insamfile  = sys.argv[1]
    outsamfile_singlemappers = sys.argv[2]
    outsamfile_multimappers = sys.argv[3]    
except:
    sys.exit("Please, give .sam in and outputfiles (single and multimappers) respectively\nOptionally, give two suffixes in order to split reads 1 and 2")
    
if len(sys.argv)>=4:
    [R1_suffix, R2_suffix] = sys.argv[4].split(';')
else:
    [R1_suffix, R2_suffix] = ['','']
    
print('-- Extracting cigar and mismatch values, separating single and multimappers --')    
    
# testing purposes
#insamfile = 'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.f2.sam'
#outsamfile_singlemappers = 'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.singlemappers.mw.sam'
#outsamfile_multimappers = 'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.multimappers.mw.sam'
# [R1_suffix, R2_suffix] = ['.R1','.R2']
f     = open(insamfile)

if R1_suffix==R2_suffix:
    splitPairs = False
    print("NOT splitting into R1 and R2 files.")
else:
    splitPairs = True
    print("Splitting into R1 and R2 files.")

# set up single mappers output
# allow for either output to same file, or output to two files
outsamfile_singlemappers_RA = outsamfile_singlemappers.replace('.sam',R1_suffix+'.sam')
outsamfile_singlemappers_RB = outsamfile_singlemappers.replace('.sam',R2_suffix+'.sam')

f_out_smap=[] # make it an array to allow for splitting of pairs
f_out_smap.append(open(outsamfile_singlemappers_RA, 'w'))
if not splitPairs:
    f_out_smap.append(f_out_smap[0])
else:
    f_out_smap.append(open(outsamfile_singlemappers_RB, 'w'))
    
#f_out_mmap = open(outsamfile_multimappers, 'w')

# Set up multimapper file output
outsamfile_multimappers_RA = outsamfile_multimappers.replace('.sam',R1_suffix+'.sam')
outsamfile_multimappers_RB = outsamfile_multimappers.replace('.sam',R2_suffix+'.sam')

f_out_mmap=[] # array to allow for splitting of pairs
f_out_mmap.append(open(outsamfile_multimappers_RA, 'w'))
if not splitPairs:
    f_out_mmap.append(f_out_mmap[0])
else:
    f_out_mmap.append(open(outsamfile_multimappers_RB, 'w'))

fd = open('debug.txt','w')

actualReadFlag = False
counter=0; aR_counter=0; nr_unmapped=0; nr_mapped=0
nr_singlemap=0;nr_multimap=0;nr_paired=0;nr_singleend=0
anomaly_flag=False;nr_anomalies=0
while (True):

    counter+=1
    CG = ['','']
    
    if anomaly_flag:
        line1=line2
        anomaly_flag = False
    else:
        line1 = f.readline()
        
    #print('R1')
    
    if (line1==''):
        # print('eof')
        break
    
    # ignore first part of file which doesn't contain reads
    if not actualReadFlag:
        if (re.match('^@', line1)):
            f_out_smap[0].write(line1)
            f_out_mmap[0].write(line1)
            if splitPairs:
                f_out_smap[1].write(line1)
                f_out_mmap[1].write(line1)
            continue
        else:
            # set signal that we'll now be reading actual reads
            actualReadFlag = True
            #print('actual reads starting')
            
    #print('==============================')
    aR_counter += 1 # actual read (pairs)        
           
    # For debug purposes
    #if counter>250: 
    #    break
    
    # now determine if these are paired reads, assuming 1st read is representative
    line1_split = line1.split('\t') # use the flag for this
    FLAG1 = line1_split[1]
    if (bin(int(FLAG1))[-1]=='1'):
        paired = True
        #print('paired-end')
    else:
        paired = False
        #print('single-end')

    # If self unmapped or mate unmapped, let's just skip this one
    if (bin(int(FLAG1))[-3]=='1') or (bin(int(FLAG1))[-4]=='1'):
        nr_unmapped += 1
        if paired:
              line2 = f.readline() # might as well skip its mate also
              
              # FOR DEBUGGING PURPOSES
              line2_split = line2.split('\t')  # DEBUGGING ONLY
              
              if not (line1_split[0] == line2_split[0]):
                  nr_anomalies += 1
                  fd.write('>>>> ANOMOLOUS unmapped read detected, claimed pair, but didnt find both'+'\n')
                  fd.write(line1_split[0]+'\n')
                  #print(line1)
                  #print(line2)
                  anomaly_flag = True
                  continue
              # DEBUGGING END
              
              #print('R2')
              #print('Unmapped, discarded:')
              #print('l1===>'+line1)
              #print('l2===>'+line2)
        continue
    else:
        nr_mapped += 1

    # Now identify which read is which
    if paired:
        if (bin(int(FLAG1))[-7]=='1'):
            i1=0; i2=1
        else:
            i1=1; i2=0
        # also keep count # pair/single reads
        nr_paired += 1
    else:
        i1=0
        # also keep count # pair/single reads
        nr_singleend += 1
            
    # Determine whether this is a singlemapper or multimapper; 
    NH = int([substr for substr in line1_split if re.match('NH:i:', substr)][0][5:])
                
    # Now extract the CIGAR string
    CG[i1] = line1_split[5]
        
    # do stuff differently depending read is paired or not
    if paired:
        
        # Read the mate (note: input file should be sorted in right way)
        line2 = f.readline()
        #print('R2')
        
        # FOR DEBUGGING PURPOSES
        line2_split = line2.split('\t')  # DEBUGGING ONLY
      
        if not (line1_split[0] == line2_split[0]):
            nr_anomalies += 1
            fd.write('>>>> ANOMOLOUS MAPPED read detected, claimed pair, but didnt find both'+'\n')
            fd.write(line1_split[0]+'\n')
            #print(line1)
            #print(line2)
            anomaly_flag = True
            continue
        # DEBUGGING END
        
        # print('Pair found')

        # search for nM (not NM!), and extract value
        # (since this is an optional field, opted for searching instead taking defined column)
        nM = [substr for substr in line1_split if re.match('nM:i:', substr)][0][5:]
        
        # Also read line and determine cigar here
        line2_split = line2.split('\t')                
        CG[i2] = line2_split[5]
        
        # print('nM: '+nM)
        # print('CG:' + CG[0] + '/'+CG[1])
    
        # Now export the two reads to output file
        # Use the i1 and i2 -- which identify which read is 1 and which one 2 --
        # to put R1 and R2 in their respective files
        if (NH==1):
            #print('writing single-mapper')
            f_out_smap[i1].write(line1_split[0]+';nM:'+nM+';C1:'+CG[0]+';C2:'+CG[1]+'\t'+'\t'.join(line1_split[1:]))
            f_out_smap[i2].write(line2_split[0]+';nM:'+nM+';C1:'+CG[0]+';C2:'+CG[1]+'\t'+'\t'.join(line2_split[1:]))   
            #print(line1_split[0])
            #print(line2_split[0])
            nr_singlemap += 1
            
        elif (NH>1):
            #print('writing multi-mapper')            
            f_out_mmap[i1].write(line1_split[0]+';nM:'+nM+';C1:'+CG[0]+';C2:'+CG[1]+'\t'+'\t'.join(line1_split[1:]))
            f_out_mmap[i2].write(line2_split[0]+';nM:'+nM+';C1:'+CG[0]+';C2:'+CG[1]+'\t'+'\t'.join(line2_split[1:]))
            #print(line1_split[0])
            #print(line2_split[0])
            nr_multimap += 1

        else:
            sys.exit('Unexpected NH value.')
            break
    
    else:
        
        NM = [substr for substr in line1_split if re.match('NM:i:', substr)][0][5:]
        # print(str(NM))
        
        # print('Single end read found.')
        
        # Now export the single read to output file
        if (NH==1):
            f_out_smap[0].write(line1_split[0]+';NM:'+NM+';CG:'+CG[0]+CG[1]+'\t'+'\t'.join(line1_split[1:])) # note: if single read, one of CG will be ''
            nr_singlemap += 1
        elif (NH>1):
            f_out_mmap[0].write(line1_split[0]+';NM:'+NM+';CG:'+CG[0]+CG[1]+'\t'+'\t'.join(line1_split[1:]))
            nr_multimap += 1
        else:
            sys.exit('Unexpected NH value.')
            break
    
    if not (line1_split[0] == line2_split[0]):
        print('major error')
        break
    # REMOVE
    HI1 = [substr for substr in line1_split if re.match('HI:i:', substr)][0][5:]
    HI2 = [substr for substr in line2_split if re.match('HI:i:', substr)][0][5:]
    if not (HI1 == HI2):
        print('major error: HI didnt match')
        break


f.close()
f_out_smap[0].close()
f_out_mmap[0].close()
if splitPairs:
    f_out_smap[1].close()
    f_out_mmap[1].close()
fd.close()

print('Done. Summary:')

print('===')

print('number of mapped read (pairs) processed: '+str(nr_mapped))
print('number of unmapped processed: '+str(nr_unmapped))

print('number of single-mappers: '+str(nr_singlemap))
print('number of multi-mappers: '+str(nr_multimap))

print('number of paired end reads: '+str(nr_paired))
print('number of single end reads: '+str(nr_singleend))

print('number of anomalies: '+str(nr_anomalies))

print('===')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    