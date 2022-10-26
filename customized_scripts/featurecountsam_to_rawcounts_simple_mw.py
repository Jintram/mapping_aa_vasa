#!/usr/bin/env python3
"""
Created on Fri Sep 30 14:11:44 2022

@author: m.wehrens
"""
import sys

# Important notes
# This scripts is based on two assumptions:
# - That in the first column of the file, the 7th element is the well name
# - That the final column of the file contains the gene annotation   

inputfile = sys.argv[1]
outputfile = sys.argv[2]

# 
# inputfile = "/Volumes/workdrive_m.wehrens_hubrecht/temporary_files/Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT.nonRibo_E99_Aligned.out.bam.featureCounts.sam"
# outputfile = "/Volumes/workdrive_m.wehrens_hubrecht/temporary_files/Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT.counts-raw.tsv"

# inputfile = "/Volumes/workdrive_m.wehrens_hubrecht/temporary_files/Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT.nonRibo_E99_Aligned.out.bam.featureCounts.sam"
# outputfile = "/Volumes/workdrive_m.wehrens_hubrecht/temporary_files/Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT.counts-raw.tsv"

cnt={}

# First process the file per line
with open(inputfile) as f:
    for i, line in enumerate(f):    
        
        if not (line[0]=='@'):
            
            # process the line
            splitline = line.rstrip().rsplit('\t')
            infoline = splitline[0] 
            
            # Select the well nr, which is #7 in the infoline
            CB   = infoline.rsplit(';')[6].rsplit(':')[1]
            # Select the gene, which is the last slot of the sam file
            gene = splitline[len(splitline)-1].rsplit(':')[2]
            
        
            # add 1 to the count in the table for every observation
            try:
                cnt[CB][gene]+=1
            except:
                try:
                    cnt[CB][gene]=1
                except:
                    cnt[CB] = {}
                    cnt[CB][gene]=1
            
        
        
print('Done creating count table, '+str(i)+' reads processed')
        

# Now output the file
with open(outputfile, 'w') as the_file:
    the_file.write('gene\tcell\tcount\n')
    for well in cnt.keys():
        #print(well)
        for gene in cnt[well]:
            #print(well+'\t'+gene+'\t'+str(cnt[well][gene]))
            the_file.write(gene+'\t'+well+'\t'+str(cnt[well][gene])+'\n')
        

print('Done, output file written, goodbye..')        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        