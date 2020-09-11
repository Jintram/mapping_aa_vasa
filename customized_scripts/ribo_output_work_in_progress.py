#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:44:39 2020

@author: m.wehrens
"""

import numpy as np

for read in samfile.fetch('chr1', 100, 120):
     print read
     
 for i, r in enumerate(bam.fetch(until_eof = True)):
     print(r)
     break
     
# let's run a little check 
bam = pysam.AlignmentFile(bamfile)
for i, r in enumerate(bam.fetch(until_eof = True)):
    #print(i)
    
    if np.mod(i,2)==0:
        r1=r

    if np.mod(i,2)==1:
        r2=r
    
        if (r1.qname==r2.qname):
            print('yas')
        else:
            print('NOOO!')
            sys.exit()
    
     
