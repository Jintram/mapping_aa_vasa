#!/usr/bin/env python3
# Reads R1.fastq and R2.fastq files;
# selects reads with proper cell barcode;
# produces a new _cbc.fastq.gz file.
import sys, os
import itertools as it
import argparse as argp
import numpy as np
import gzip
import pandas as pd
from pandas.io.parsers import read_csv
from collections import Counter
import glob
import re

#### function to identify cells from barcodes, allowing some edit distances ####
def find_compatible_barcodes(barcode, HDmax = 0):
    """Given a barcode sequence and a maximum Hammin distance, it returns a list of compatible barcode sequences"""
    nt = ['N'] if HDmax == 0 else ['N','C','T','G','A']
    HDmax = 1 if HDmax == 0 else HDmax

    compatible_barcodes = set([barcode])
    for hd in range(1, HDmax+1):
        comb = [''.join(l) for l in it.product(nt, repeat = hd)]
        for c in comb:
            for p in it.permutations(range(len(barcode)), hd):
                s0 = barcode
                for x, l in zip(p, c):
                    s0 = s0[:x] + l + s0[x+1:]
                compatible_barcodes.add(s0)
    return list(compatible_barcodes)

#### check input variables ####
parser = argp.ArgumentParser(description = 'Concatenates bcread to bioread qname.')
parser.add_argument('--fqf', help = 'Fastq files names, without _Rx.fastq.gz')
parser.add_argument('--bcread', '-bcr', help = 'read where to find the barcode (umi+cell)', choices = ['R1', 'R2'], default = 'R1')
parser.add_argument('--bioread', '-bior', help = 'read where to find biological information', choices = ['R1', 'R2'], default = 'R2')
parser.add_argument('--lencbc', '-lcbc', help = 'cell barcode length (integer)', type = int, default = 8)
parser.add_argument('--lenumi', '-lumi', help = 'umi length (integer)', type = int, default = 6)
parser.add_argument('--umifirst', help = 'logical variable: umi before cel barcode', action = 'store_true')
parser.add_argument('--cbcfile', '-cbcf', help = 'cell specific barcode file. Please, provide full name')
parser.add_argument('--cbchd', help = 'collapse cell barcodes with the given hamming distance', type = int, default = 0)
parser.add_argument('--outdir', help = 'output directory for cbc.fastq.gz and log files', type = str, default = './')
args = parser.parse_args()

fqr = args.fqf
bcread = args.bcread
bioread = args.bioread
lcbc = args.lencbc
lumi = args.lenumi
umifirst = args.umifirst
cbcfile = args.cbcfile
hd = args.cbchd
outdir = args.outdir
target_primer_len = 20

# wdir = os.getcwd() # allows for easier interaction during testing

#### Find input fastq files ####
fq1 = glob.glob(fqr + '_R1*.fastq.gz')
fq2 = glob.glob(fqr + '_R2*.fastq.gz')
print(fq1, fq2)

if len(fq1) != 1 or len(fq1) != 1:
    sys.exit("Please, be more specific with fqf tag")

fq1 = fq1[0]; fq2 = fq2[0]

if not os.path.isfile(fq1) or not os.path.isfile(fq2):
    sys.exit('fastq files not found')

#### Read barcodes and expand set according to input hamming distance ####
if not os.path.isfile(cbcfile):
    sys.exit("Barcode file not found")

bc_df = read_csv(cbcfile, sep = '\t', names = ['bc','cellID'], index_col = 0)
bc_df['compatible_bcs'] = bc_df.apply(lambda x: find_compatible_barcodes(x.name, hd), axis = 1)
cnt_allbcs = Counter([x for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs']])
allbc_df = pd.DataFrame({x: {'cellID': bc_df.loc[idx,'cellID'], 'original': idx} for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs'] if cnt_allbcs[x]==1}).T

### Create output directory if it does not exist ####
if not os.path.isdir(outdir):
    os.system('mkdir '+outdir)

#### Read fastq files and assign cell barcode and UMI ####
fout_R1 = open(outdir + '/' + fqr + '_R1_cbc.fastq', 'w')
fout_R2 = open(outdir + '/' + fqr + '_R2_cbc.fastq', 'w')

ns = 0
readcount = 0
reads_polyT = 0
reads_targeted = 0
reads_unclassified = 0
with gzip.open(fq1) as f1, gzip.open(fq2) as f2: 
    for idx, (l1, l2) in enumerate(zip(f1, f2)):

        # read current line from R1, R2
        l1, l2 = str(l1.rstrip().rsplit()[0], 'utf-8'), str(l2.rstrip().rsplit()[0], 'utf-8')
        l = np.mod(idx,4)

        # Each entry consists of 4 consecutive lines, go over them 1-by-1
        # line 1, read to nX ("name"), check consistency
        if l == 0:
            n1, n2 = l1, l2
            if not n1 == n2:
                print (n1, n2)
                sys.exit('fastq files not synchronized (@name)')
        # line 2, read to sX ("sequence")
        if l == 1:
            s1, s2 = l1, l2
        # line 3, read to pX ("plus sign", separator)
        if l == 2:
            p1, p2 = l1[0], l2[0]
            if not p1 == p2: # == '+':
                print(l1, l2, p1, p2)
                sys.exit('fastq files not synchronized (+)')
        # line 4, read to qX ("quality"), also start processing
        if l == 3:
            readcount += 1
            classification = "none"
            extra_trimming = 0
        
            q1, q2 = l1, l2
            if len(q1) != len(s1) or len(q2) != len(s2):
                sys.exit('phred and read length not mathch!')

            # all 4 lines read, now we can process

            # customized pre-processing of the read, to look for the targeted
            # primer
            
            # if re.search("|".join(target_primers), s1): # not robust against artefacts
            # determine position primer in read
                #target_primer_pos = re.search("|".join(target_primers), s1).span()
                #extra_trimming = target_primer_pos[1]-target_primer_pos[0]
                
            # get the sequence at the position where targeted primer should be located if present
            # (otherwise it'll be poly-T)
            # Note that this cannot deal with primers of different lengths yet,
            # this can of course be adjusted, probably easiest by writing
            # a little function that goes over each primer
            # For now we won't do that, since that'll cost more processing time
            seq_at_primer_pos = s1[lumi+lcbc:lumi+lcbc+target_primer_len]
            # now act accordingly
            if (seq_at_primer_pos in target_primers):     
                # assign primer sequence as classification
                classification = seq_at_primer_pos
                # set extra trimming to remove the primer sequence
                extra_trimming = target_primer_len 
                # message
                print("Classified as targeted read ("+classification+")")
                reads_targeted += 1
            # classify as poly-T when 20 Ts found (primer has 26, but not sure reliable poly-read); 26=TTTTTTTTTTTTTTTTTTTTTTTTTT
            # poly-T sequences will be removed later by trim galore / cutadapt
            elif (seq_at_primer_pos[0:10] == 'TTTTTTTTTT'):                 
                classification = "polyT"
                reads_polyT += 1
                print("Classified as poly-T read.")                
            else:
                reads_unclassified += 1
                print("Read not classified")
                print(s1)
                
            # determine BC + UMI
            if bcread == 'R1':
                bcseq = s1[:lumi+lcbc]
                bcphred = q1[:lumi+lcbc]
                s1 = s1[lumi+lcbc+extra_trimming:]
                q1 = q1[lumi+lcbc+extra_trimming:]
            elif bcread == 'R2':
                bcseq = s2[:lumi+lcbc]
                bcphred = q2[:lumi+lcbc]
                s2 = s2[lumi+lcbc+extra_trimming:]
                q2 = q2[lumi+lcbc+extra_trimming:]  
                
            if not umifirst:
                cellbcseq = bcseq[:lcbc]
                umiseq = bcseq[lcbc:]
                cellbcphred = bcphred[:lcbc] # quality string
                umiphred = bcphred[lcbc:]  # quality string
            else:
                cellbcseq = bcseq[lumi:]
                umiseq = bcseq[:lumi]
                cellbcphred = bcphred[lumi:]  # quality string
                umiphred = bcphred[:lumi]  # quality string

            try:
                # get cell ID from BC (failures are skipped I guess)
                cellID, originalBC = allbc_df.loc[cellbcseq]
                ns += 1
                cellbcphred = ''.join([chr(ord(c)+32) for c in cellbcphred])                
                umiphred = ''.join([chr(ord(c)+32) for c in umiphred])

                # prepare new name, and output respective reads to files
                name = ';'.join([n1] + [':'.join(x) for x in zip(['SS','CB','QT','RX','RQ','SM','CL'], [cellbcseq, originalBC, cellbcphred, umiseq, umiphred, str(cellID).zfill(3), classification])])
                s, q = (s1, q1) if bioread == 'R1' else (s2, q2)
                     
                fout_R1.write( '\n'.join([name, s1, '+', q1, '']))
                fout_R2.write( '\n'.join([name, s, '+', q, '']))
                
            except: 
                print('Fail - skip')
                continue
            
        if readcount > 1000:
            break

nt = (idx+1)/4
fout_R1.close()
fout_R2.close()

#### LOG ####
fout = open(outdir + '/' + fqr + '.log', 'w')
fout.write('=> to generate cbc file <=\n')
fout.write(', '.join(['fastq file:', str(fqr),'\n']))
fout.write(', '.join(['full barcode in:', str(bcread),'\n']))
fout.write(', '.join(['biological read in:', str(bioread), '\n']))
fout.write(', '.join(['cell specific barcode length:', str(lcbc), '\n']))
fout.write(', '.join(['umi length:', str(lumi), '\n']))
fout.write(', '.join(['umi goes first:', str(umifirst),'\n']))
fout.write(', '.join(['total sequenced reads:', str(nt), '\n']))
fout.write(', '.join(['reads classified as poly-T:', str(reads_polyT), '\n']))
fout.write(', '.join(['reads classified as targeted:', str(reads_targeted), '\n']))
fout.write(', '.join(['reads not classified:', str(reads_unclassified), '\n']))
fout.write(', '.join(['reads with proper barcodes:', str(ns), str(1.0*ns/nt), '\n']))

fout.close()

#### zip fastq file ####
os.system('gzip '+ outdir + '/' + fqr + '_cbc.fastq')
