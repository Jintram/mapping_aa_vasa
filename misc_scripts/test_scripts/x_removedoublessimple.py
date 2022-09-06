#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 22:07:08 2020

@author: m.wehrens
"""

import sys

# Just removing double lines


infile  = 'TS_MYBPC_primer_reads_by_grep.txt'
outfile = 'TS_MYBPC_primer_reads_by_grep_nodouble_simple.txt'

fi = open(infile)
fo = open(outfile, 'w')

pair1_A = fi.readline()
pair1_B = fi.readline()

previous_name_A = pair1_A.split(sep='\t')[0]
previous_name_B = pair1_B.split(sep='\t')[0]

if (previous_name_A != previous_name_B):
        print('names of assumed pair do not match')
        sys.exit()

fo.write(pair1_A)
fo.write(pair1_B)


while(True):

    # read new pair
    pair1_A = fi.readline()
    # check for end of file
    if (pair1_A == ''):
        break
    pair1_B = fi.readline()
    
    # determine read name
    current_name_A = pair1_A.split(sep='\t')[0]
    current_name_B = pair1_B.split(sep='\t')[0]    
    
    # sanity check
    if (current_name_A != current_name_B):
        print('names of assumed pair do not match')
        print(current_name_A)
        print(current_name_B)
        sys.exit()
    
    if (previous_name_A != current_name_A):
        fo.write(pair1_A)
        fo.write(pair1_B)
    
    previous_name_A = current_name_A
    
fi.close()
fo.close()










