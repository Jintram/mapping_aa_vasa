#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:30:26 2020

@author: m.wehrens
"""

path = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/testrun3/extracting_seq_test/'

bedsingle = path+'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.R2.nsorted.singlemappers_genes.bed'
bedmulti = path+'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.R2.nsorted.multimappers_genes.bed'

bedsingle_R1 = 'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.R1.nsorted.singlemappers_genes.bed'
bedmulti_R1 = 'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.R1.nsorted.multimappers_genes.bed'

output = path+'S1_test2_'
protocol = 'vasa'