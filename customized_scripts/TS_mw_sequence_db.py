#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 16:11:57 2020

@author: m.wehrens
"""

###############################################################################

import pysam

# lib='MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat'

# The bamfile contains the reads with their mapping
# (presented as separate entries for the two reads, 
# and also separate for multi-mappers)
bamfile = lib+'_TS.nonRibo_E99_Aligned.nsorted.out.bam'
# The decisionfile contains information about to which
# gene this read was assigned by Anna's rules
decisionfile = lib+'_merged_decisions_sorted.tsv'

###############################################################################

# open the bam file
bam = pysam.AlignmentFile(bamfile)
# open the decision file
f=open(decisionfile)
    
# similar strategy as anna, 
# Loop over the reads, collect them,
# but here the idea is that we'll then 
# use the name to look up to which gene it was finally mapped
# using the decisions file, and couply those two quantities
# we aim to create an output file with each line:
# gene_mapped_to|read1or2|sequence|barcode|primer|qname

for i, bam_entry in enumerate(bam.fetch(until_eof = True)):
    
    print(i)
    
    # init in first loop
    if i==0:
        
        previous_name = bam_entry.qname
        current_read_collection = [bam_entry]
    
    else:
        
        # if this mapping still belongs to this read;
        # add to collection
        if bam_entry.qname==previous_name:
            current_read_collection.append(bam_entry)
            
        # if not, start processing the read
        else:
            decision_entry = f.readlines(1)
            decision_rname = decision_entry[0].split('\t')[0]
            
            current_read_collection[0].qname == decision_rname
            
            break
    
    #print(bam_entry.qname)

    
    if i==3: break




f.close()













