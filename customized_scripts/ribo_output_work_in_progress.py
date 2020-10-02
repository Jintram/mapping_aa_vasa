#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:44:39 2020

@author: m.wehrens
"""


# A little note:
# 
# One can find mates for indexed bam files by
# calling:
# bam.mate(r1) 
# where "bam" is created with:
# bam = pysam.AlignmentFile(bamfile)
# 
# The file we're using is however not indexed, but 
# sorted by name.
# For the paired stuff, it contains:
# r1a, r1b, r2a, r2b consecutively,
# where I assume a and b respectively come from different 
# mapping


import numpy as np

for read in samfile.fetch('chr1', 100, 120):
     print(read)
     
 for i, r in enumerate(bam.fetch(until_eof = True)):
     print(r)
     break
     
# let's run a little check 
bam = pysam.AlignmentFile(bamfile)
for i, r in enumerate(bam.fetch(until_eof = True)):
    print(i)
    
    if np.mod(i,4)==0:
        r1=r

    if np.mod(i,4)==1:
        r2=r
    
    if np.mod(i,4)==2:
        r3=r

    if np.mod(i,4)==3:
        r4=r

    
        if (r1.qname==r2.qname):
            print('r1=r2')
        else:
            print('NOOO!')
            sys.exit()
            
        if (r1.qname==r2.qname==r3.qname==r4.qname):
            print('All 4 are same')
        else:
            print('no 4 same!')
            sys.exit()
            
            
            
            
            
            
            
            
    
     
