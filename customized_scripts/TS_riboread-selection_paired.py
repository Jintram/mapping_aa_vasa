#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
import pysam
from collections import Counter

# check input
try:
    bamfile = sys.argv[1]
    stranded = sys.argv[2]
    output = sys.argv[3]
except:
    sys.exit("Please, give: \n(1) input bam file;\n(2) stranded (y/n);\n(3) prefix for output file")

bam = pysam.AlignmentFile(bamfile)
bam_out = pysam.AlignmentFile(output + '.Ribo.bam', template = bam, mode = 'wb')
nreads = 0
nunmapped = 0
nmapped = Counter()
nmapped_ = 0
fout_R1 = open(output + '_R1.nonRibo.fastq', 'w')
fout_R2 = open(output + '_R2.nonRibo.fastq', 'w')
reads1 = []
reads2 = []
# go over bam file entries
# note that each entry is duplicated because it was mapped twice
for i, r in enumerate(bam.fetch(until_eof = True)):

    # initialize "reads" during first loop    
    if i == 0:

        # note that "reads" will be re-initialized every time 
        # a series of entries of same read is broken (see below)
        reads = [r]
        nreads += 1

    else:
        # while the name of the read is the same,
        # add it to the collection of reads
        if r.qname == reads[0].qname:
            reads.append(r)
            
        else: 
            # else, start a new series (see below), but first process this series
            # of reads
            
            # if all are unmapped, this is not a ribosome read,
            # and a single entry of the read can be written to
            # the output fastq file
            if all([rs.is_unmapped for rs in reads]): # all umapped => fastq file
                
                # output to fastq file
                # ===
                
                # in theory, bam.mate(r0) would give the mate, however, 
                # -- see https://github.com/pysam-developers/pysam/issues/472 --
                # due to an issue, this won't work if they're unmapped because
                # r.rname and r.mrnm (the names of the ref seq of the read and mate's read)
                # are not set, which happens if they're not mapped.
                # 
                # .. but luckily we can do this manually
                
                # determine R1
                r1=[r for r in reads if r.is_read1][0]
                # determine R2
                r2=[r for r in reads if r.is_read2][0]
                
                # first pre-proces read for output
                r1.qname = '@' + r1.qname
                r2.qname = '@' + r2.qname

                # admin
                nunmapped += 1                
                                
                # write to file in fastq format
                fout_R1.write('\n'.join([r1.qname, r1.get_forward_sequence(), '+', r1.qual]) + '\n')
                fout_R2.write('\n'.join([r2.qname, r2.get_forward_sequence(), '+', r2.qual]) + '\n')
                    # get_forward_sequence retrieves original sequence (otherwise 
                    # sequences are sometimes reverse complements)
                                        
                # original code by Anna for handling non-paired reads
                #r0.qname = '@' + r0.qname
                #fout.write('\n'.join([r0.qname, r0.seq, '+', r0.qual]) + '\n')
                #nunmapped += 1
                
            else:
                
                if stranded == 'n': 
                    mapreads = [x for x in reads if not x.is_unmapped]
                    
                # ignore reverse hits
                elif stranded == 'y':
                    mapreads = [x for x in reads if (not x.is_unmapped) and (not x.is_reverse)]
                    
                if len(mapreads) >= 1: # at least one is mapped properly => bam file for ribo data
                    rgtag = '_'.join(sorted([r.get_tag('RG').rsplit('.')[-1] if 'RG' in [t[0] for t in r.get_tags()] else '-' for r in mapreads]))
                    rgtag = rgtag.replace('-ribo','')
                    nmapped.update([rgtag])
                    nmapped_ += 1 # redundant with nmapped
                    r0 = mapreads[0]
                    new_tags = [t if t[0] != 'RG' else ('RG',rgtag) for t in r0.tags]
                    r0.tags = new_tags
                    bam_out.write(r0)
                    
                else:
                    # if after selection (ie ignoring reverse) no reads are mapped 
                    # it's still a non-ribo read, and we can output it.
                    
                    # output to fastq file
                    # ===                                        
                    
                    # determine R1
                    r1=[r for r in reads if r.is_read1][0]
                    # determine R2
                    r2=[r for r in reads if r.is_read2][0]
                    
                    # first pre-proces read for output
                    r1.qname = '@' + r1.qname
                    r2.qname = '@' + r2.qname
    
                    # admin
                    nunmapped += 1                
                                    
                    # write to file in fastq format
                    fout_R1.write('\n'.join([r1.qname, r1.get_forward_sequence(), '+', r1.qual]) + '\n')
                    fout_R2.write('\n'.join([r2.qname, r2.get_forward_sequence(), '+', r2.qual]) + '\n')
                        # get_forward_sequence retrieves original sequence (otherwise 
                        # sequences are sometimes reverse complements)
                                           
                    # original code by Anna
                    #r0 = reads[0]
                    #r0.qname = '@' + r0.qname
                    #fout.write('\n'.join([r0.qname, r0.seq, '+', r0.qual]) + '\n')
                    #nunmapped += 1
                    
            # since the new read's name didn't match with the current
            # series of reads, start a new series
            reads = [r]
            nreads += 1

fout_R1.close()
fout_R2.close()
bam.close()
bam_out.close()

fout = open(output + '.ribo-map.log', 'w')
fout.write('Number of reads: ' + str(nreads)+'\n')
fout.write('Number of unmapped reads: '+str(nunmapped)+'\n')
fout.write('Number of mapped reads: '+'\n')
for tag in nmapped:
    fout.write('\t'+tag+': '+str(nmapped[tag])+'\n')
fout.close()
os.system('gzip '+output+'.nonRibo.fastq')
