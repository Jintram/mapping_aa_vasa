#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 16:19:43 2020

@author: m.wehrens
"""

# this script will now merge all cells into one dataframe

import pandas as pd
from functools import reduce
import numpy as np

f = open('percell/cells_with_reads.tsv')

df_list = []
colnames = []
for line in f:

    df_list.append(pd.read_table('percell/rsem_'+line.split('\n')[0]+'.genes.results', sep='\t', index_col = 0))
    colnames.append(line.split('\n')[0])
    
f.close()
A_test=df_list[0]
    
# Create list of data frames that will later serve as columns of output
df_list_sub = [df.iloc[:,[3]] for df in df_list]

# Rename the columns to respective samples
for idx in range(len(df_list_sub)):
    df_list_sub[idx].rename(columns={'expected_count': colnames[idx]}, inplace=True)

# Now merge the dataframes
df_merged = df_list_sub[0].join(df_list_sub[1:], how='outer')

# Now export the merged dataframe as tsv file
df_merged.to_csv('RSEM_countTable.csv', sep='\t', na_rep='')


# HOWEVER NOW WE HAVE RAW READCOUNTS, BUT NOT TAKEN INTO ACCOUNT UMIs...
# Oops..









