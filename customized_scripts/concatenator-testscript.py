#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 15:15:21 2020

@author: m.wehrens
"""

fqr = 'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001'
bcread = 'R1'
bioread = 'R2'
lcbc = 8
lumi = 6
umifirst = True
cbcfile = '/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/bc_celseq2.tsv'
hd = 0
outdir = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3'

#target_primers = ['GGTCATGAGTCCTTCCACGA','AACGTACAAAGTGGGGATGG','ACCCCTGGAGACTTTGTCTC',
#                  'TGCGCTCCAGGATGTAGCCC','CCATCGTAGGCAGGCGGCTC','CAGGATGTAGCCCAGGATGG']

primerfile = '/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/targetedprimers.tsv'

#os.chdir('/Volumes/workdrive_m.wehrens_hubrecht/mapping-test/pilot3')
#wdir = '/Volumes/workdrive_m.wehrens_hubrecht/mapping-test/pilot3/'

fq1 = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001_R1_001.fastq.gz'
fq2 = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001_R2_001.fastq.gz'

