#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 16:59:12 2020

@author: m.wehrens
"""

# To run count table file locally for testing purposes:
bedsingle = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.singlemappers_genes.bed'
bedmulti = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.multimappers_genes.bed'
output = 'test_count_out'
protocol = 'smartseq_UMI'

# for single read test
bedsingle = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001_SS_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.singlemappers_genes.bed'
bedmulti = '/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001_SS_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.multimappers_genes.bed'
