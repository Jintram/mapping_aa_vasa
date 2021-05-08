#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 16:44:34 2021

@author: m.wehrens
"""

# This is a very simple conversion script, converting from .bed
# format to Anna's custom format.
#
# Final columns should be:
# chromosome_nr start_pos end_pos strand ENSEMBLID_GENENAME_GENEBIOTYPE_EXONORINTRON GENE_LENGTH GENE_START GENE_END

path_in = '/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/downloaded_ensembl/81/head_Homo_sapiens.GRCh38.81.bed'
path_out = '/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/downloaded_ensembl/81/custom_head_Homo_sapiens.GRCh38.81.bed'

# dict mapping the columns
the_fields = ['chromosome_nr','start_pos','end_pos','','','','','','',]
name2idx = {'chromosome_nr':1, 
            'start_pos':2,
            'end_pos':3,
            'strand':6,            
            'ensemblid':4,
            'genename':21,
            'biotype':25,
            'exonorgene':8
            # gene length should be calculated (depends strand!)
            # gene start and gene end are infered from 'gene' annotation rows
            }
name2idx.keys()
# 

f_in  = open(path_in,'r')
f_out = open(path_out,'w')

while (True):
    
    line = f_in.readline()
    
    columns = line.split()
    
    line_formatted = {'chromosome_nr':columns[name2idx['chromosome_nr']]}
    
    {'chromosome_nr':columns[name2idx['chromosome_nr']] for X in the_fields}

    if current_group.append()
    
    
    
    
    
    



close(f_in)
close(f_out)