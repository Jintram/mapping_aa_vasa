#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 11:37:29 2021

This is a simple file that loops over a sam file, and splits the file into multiple
parts according to the cell bar codes.

Probably there's faster ways to do this, but for now let's try this.
(Now took ±an hour @ HPC for 60 GB bam file / 260 GB sam file.)

@author: m.wehrens
"""

#base_folder='/Volumes/workdrive_m.wehrens_hubrecht/mapping-test/wang_n2/'
base_folder = '/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/mapping_may25_download/'
input_file = 'GSM2970358_N2_LV_cat_nc.nonRibo_E99_Aligned.out'

from pandas.io.parsers import read_csv
import numpy as np
import time


# First import the barcodes.
#cbcfile = '/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/bc_takara_map_tcr_v1.tsv'
cbcfile = '/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04/bc_takara_map_tcr_v1.tsv'
bc_df = read_csv(cbcfile, sep = '\t', names = ['bc','cellID'], index_col = 0)
    # len(bc_df.index.values)
    
# Then create a simple dictionary to map respective reads to split the barcodes 
# over files
all_barcodes = bc_df.index.values
part_dict_ = np.repeat(range(10),np.ceil(len(all_barcodes)/10))
the_assignment = {all_barcodes[idx]:part_dict_[idx] for idx in range(0,len(all_barcodes))}
    # len(part_dict_)
    # len(all_barcodes)

# Time this operation
start = time.time()

# Now open all 10 files
files = [open(base_folder+'part.'+str(idx)+'.'+input_file+'.sam', mode='w') for idx in range(10)]

# Open input file
inputFile = base_folder+input_file+'.sam'
f_input = open(inputFile, mode='r')

# Now read the input file
"""
for line in f_input:

    # identify bc
    current_BC = line.split(';CB:')[1].split(';QT:')[0]
    
    # write line to according file
    files[the_assignment[current_BC]].write(line)
"""

total_line_count=0
while f_input:
    
    line  = f_input.readline()
    
    # header lines (for speedup could be moved to separate pre loop)
    while line[0]=='@':
        
        print('writing header lines ..')
        
        for idxf in range(10):
            files[idxf].write(line)
        total_line_count += 1            
        
        line  = f_input.readline()
    
    for i in range(1000000):

        try:
            # identify bc
            current_BC = line.split(';CB:')[1].split(';QT:')[0]

            # write line to according file
            files[the_assignment[current_BC]].write(line)
            total_line_count += 1

        except:
            print('Failed to process the following line:')
            print('========')
            print(line)
            print('========')            
            raise ValueError('Read error.')
        
        
        # read line
        line  = f_input.readline()
        
        if line == "":
            break
    
    print(str(total_line_count)+' lines processed')    

    if line == "":
        break


# Now close all files
f_input.close()
for f in files:
    f.close()

# Do some final stuff
end = time.time()
print('Seconds elapsed:')
print(end - start)

print('Total lines processed: '+str(total_line_count))

print('Don\'t forget to convert to .bam')

# Took 0.36 seconds for 30.7 mb
# --> for 57 GB = 57000 mb, 57000/31=1838, 1838 times as long, or 11 minutes
# --> should be OK
























