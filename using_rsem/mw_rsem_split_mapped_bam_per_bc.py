#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 17:41:40 2020

@author: m.wehrens
"""

import re, os

insamfile = 'star_MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_Aligned.toTranscriptome.out.f2.sam'
outdir    = 'percell/'

if not os.path.isdir(outdir):
    os.mkdir(outdir)

# Strategy:

f = open(insamfile)

# The "SM" field signals to which cell it has been mapped, "CB"
# the detected bar code.

# first read through header
is_header=True
headerlines=[]
while(is_header):
    line = f.readline()
    if not (re.match('^@',line)):
        is_header=False
    else:
        headerlines.append(line)

# now sort out reads
line_part1_split = line.split('\t')[0].split(';')
# determine where field of interests are
SM_index = [i for i in range(1,len(line_part1_split)) if ('SM:' in line_part1_split[i])][0] # start at 1 to exclude name
CB_index = [i for i in range(1,len(line_part1_split)) if ('CB:' in line_part1_split[i])][0]

# now repeat reading lines
openfilelist = []
out_f = {}
while(True):

    # split line ..
    line_part1_split = line.split('\t')[0].split(';')
    
    # .. such that we can determine cellindex and barcode
    CellIndex = line_part1_split[SM_index][3:]
    BarCode   = line_part1_split[CB_index][3:]
    
    # determine current file name
    currentcellname = CellIndex+'_'+BarCode
    
    # keep track of file names in list, 
    # open new file if necessary
    if not currentcellname in openfilelist:
        out_f[currentcellname] = open(outdir + currentcellname +'.sam','w')
        out_f[currentcellname].write(''.join(headerlines))
        
        openfilelist.append(currentcellname)
    
    # write to file
    out_f[currentcellname].write(line)    
    
    # read next line
    line = f.readline()
    
    # stop if eof
    if line=='':
        break

# a file with cells that have detected reads
f_summary = open(outdir+'cells_with_reads.tsv','w')


# close all files
for key in out_f.keys(): 
    # write to summary file
    f_summary.write(key+'\n')
    # close files
    out_f[key].close()

# close other files
f.close()
f_summary.close()



# then create 384 (or less) files like "001.CGTCTAAT.sam"
# loop over each line and put them in according files
# Find out what to do with header though..? (probably non-consequential)

# write also txt file with a list of those filenames
# to be used later to put into rsem.
        
        


        
        
        
        