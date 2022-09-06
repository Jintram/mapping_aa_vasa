#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 22:42:36 2020

@author: m.wehrens
"""

# Just for reference, quality string is one of (increasing quality):
# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

filein = 'MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS_R1_cbc_val_1_HATCG.fq'
filein = 'MW-TS-S1A-Hybr-Bst_HTMH2BGXF_S25_cat_TS_R1_cbc_val_1_HATCG.fq'
filein = 'MW-TS-S2A-TargOnly-Bst_HTMH2BGXF_S26_cat_TS_R1_cbc_val_1_HATCG.fq'
#filein = 'MW-TS-S4A-Vasa-Bst_HTMH2BGXF_S27_cat_TS_R1_cbc_val_1_HATCG.fq'

mybpc3_str = 'aggatgggctgcccgccatcgtaggcaggcggctcccac'.upper()

f = open(filein)

primer_reads_counter = 0
mybpc3_reads_counter = 0
mybpc3_longr_counter = 0

mutated_position_counter = {'A':0, 'C':0, 'T':0, 'G':0}

qual_stats = {}
qual_stats_taken = {}
for c in '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~':
    qual_stats[c] = 0
    qual_stats_taken[c] = 0

UMI_stats_mutant = dict()

loopcounter = 0
while(True):
    
    l1 = f.readline().split('\n')[0] 
    if (l1==''): # end of file
        print('eof')
        break 
    l2 = f.readline().split('\n')[0]
    l3 = f.readline().split('\n')[0]
    l4 = f.readline().split('\n')[0]
    
    l1_split = l1.split(sep=';')
    if (l1_split[-1] == 'p_seq:TGCGCTCCAGGATGTAGCCC'): # put \n here if not removed earlier! (it is now)
        #print('we have a mybpc3 primer hit')
        primer_reads_counter += 1
        
        UMI = l1_split[4][3:]
        
        # now look for expected MYBPC3 sequence
        if (l2.find(mybpc3_str[0:20]) != -1):
            mybpc3_reads_counter += 1
            
            
            search2_out = l2.find(mybpc3_str)
            if (search2_out != -1):
                #print('Longer confirmation!')
                mybpc3_longr_counter += 1
                
                mybpc3_read_sel = l2[search2_out:]
                qual_read_sel  = l4[search2_out:]
                
                if len(mybpc3_read_sel)>len(mybpc3_str):
                    #print(mybpc3_read_sel[len(mybpc3_str)-3:len(mybpc3_str)+1])
                    
                    # keep track of the base at the potentially mutated position
                    base_of_interest  = mybpc3_read_sel[len(mybpc3_str)]
                    qual_of_interest = qual_read_sel[len(mybpc3_str)]
                    
                    qual_stats[qual_of_interest] += 1
                    
                    # quality check
                    if ord(qual_of_interest) <= 65:
                        continue
                    
                    qual_stats_taken[qual_of_interest] += 1                    
                    
                    if base_of_interest in ['A','C', 'T', 'G']:
                        mutated_position_counter[base_of_interest] += 1                    
                    
                        # also keep track of UMIs for detected mutated reads
                        if (mybpc3_read_sel[len(mybpc3_str)]) == 'C':
                            
                            if not UMI in UMI_stats_mutant:
                                UMI_stats_mutant[UMI] = 1
                            else:
                                UMI_stats_mutant[UMI] += 1
                        
f.close()
        
print('mybc3 primer reads: '+str(primer_reads_counter))
print('mybc3 reads:        '+str(mybpc3_reads_counter))
print('mybc3 reads ++:     '+str(mybpc3_longr_counter))

print('mybc3 WT reads:      '+str(mutated_position_counter['T']))
print('mybc3 mutant reads:  '+str(mutated_position_counter['C']))
print('mybc3 unknown reads: '+str(mutated_position_counter['G']+mutated_position_counter['A']))










    