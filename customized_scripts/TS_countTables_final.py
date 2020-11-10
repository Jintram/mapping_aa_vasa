#!/usr/bin/env python3
import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
from collections import Counter
import pickle

try:
    bedsingle = sys.argv[1]
    bedmulti = sys.argv[2]
    output = sys.argv[3]
    protocol = sys.argv[4]
except:
    sys.exit("Please, provide input bed with single mappers (1) bed with multimappers (2) output file (3) and protocol (4)\nOptionally, provide (5) single mappers R1 and (6) multi mappers R2")
    
if len(sys.argv)>4:
    bedsingle_R1 = sys.argv[5] # bedsingle contains R2
    bedmulti_R1 = sys.argv[6] # bedmulti contains R2
    mate_file_given = True
else:
    mate_file_given = False

###############################################################################
# Functions
    
#
def addReadCount(cnt, cell, gene, umi, label):
    try:
        cnt[cell][gene][umi].update([label])
    except:
        try:
            cnt[cell][gene][umi] =  Counter([label])
        except:
            try:
                cnt[cell][gene] = {umi: Counter([label])}
            except:
                cnt[cell] = {gene: {umi: Counter([label])}}
    return cnt

def get_cell_UMI(name, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';')[1:]}
        cell = name_info['SM']
        umi = name_info['RX']
    elif protocol == 'ramda': 
        cell =  name.rsplit('.')[0]
        umi = 'A'
    elif protocol == 'smartseq_noUMI': 
        name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';')[1:]}
        cell = name_info['SM']
        umi = 'A'
    return cell, umi

def gene_assignment_single(genes, labels, infos, covs, tlens):
    # Note that jS is a customized flag that signals how the mapping 
    # relates to intron/exon/gene boundaries [I guess --> check Anna -- mw]
    
    # if mapped to only one gene, it's simple
    if len(set(genes)) == 1 and len(set(labels)) == 1:
        
        gene = genes[0]; label = labels[0]; df = pd.DataFrame()
        
    # idem, chose intron label if >1 labels though
    elif len(set(genes)) == 1 and len(set(labels)) > 1:
        
        gene = genes[0]; label = 'intron'; df = pd.DataFrame()
        
    # now if mapped >1 gene, then .. 
    else:
        
        # sys.exit() # debug purposes
        
        # first acquire # mismatching bases
        nMs_string = [c.rsplit(';nM:')[1].rsplit(';jS:')[0] for c in infos]
        nMs_sum = [sum([int(i) for i in s.split('/')]) for s in nMs_string]
            # in case of pairs, the entries look like "1/3", with the slash
            # separating respective reads.
            # this approach with sum can handle both single and paired reads.
        
        # convert different mappings to dataframes
        df = pd.DataFrame({'genes': genes, 
            'labels': labels,
            'cigars': [c.rsplit(';nM:')[0].rsplit('CG:')[1] for c in infos],
            'nMs': nMs_sum, # [int(c.rsplit(';nM:')[1].rsplit(';jS:')[0]) for c in infos],
            'jSs':  [c.rsplit(';jS:')[-1] if ';jS:' in c else 'unknown' for c in infos],
            'covs': [float(cov) for cov in covs],
            'tlens': [int(tl) for tl in tlens]
            })
    
        # I think the idea below is to create a gene dataframe, containing
        # the mappings, from which mapping will be removed according to
        # "the rules" for which mapping is most "correct"
    
        # select subset of genes that have minimal number of mismatches
        df = df[df['nMs']==df['nMs'].min()]
        gdf = {g: df_g for g, df_g in df.groupby('genes') if '_TEC' not in g}
        
        # merge entries of same genes with same information
        for g in gdf:
            
            if len(gdf[g]) > 1 and len(set(gdf[g]['labels'])) == 1 and len(set(gdf[g]['jSs'])) == 1:
                gdf[g] = gdf[g].head(1)
                
        # if multiple labels remove intron reads [I guess]
        if len(gdf) > 1 and all([gdf[g].shape[0]==1 for g in gdf]):
            if len(set(df['labels'])) > 1: 
                fdf = df[df['labels']!='intron'] 
                gdf = {g: df_g for g, df_g in fdf.groupby('genes') if '_TEC' not in g}
                
        xg = []
        
        #
        # note that 'N' stands for alignment gap in cigar
        for g in gdf:
#            if (gdf[g].shape[0] > 1 and 'N' not in gdf[g]['cigars'].iloc[0]) or (gdf[g].shape[0] == 1 and 'N' in gdf[g]['cigars'].iloc[0]):
            if  (gdf[g].shape[0] == 1 and 'N' in gdf[g]['cigars'].iloc[0]):
                xg.append(g)
                
        for g in xg:
            del gdf[g]
            
        if len(gdf) > 1: 
            setjSs = set(pd.concat([gdf[g] for g in gdf])['jSs']); xg = []
            if len(setjSs) > 1 and 'IN' in setjSs: 
                xg = [g for g in gdf if 'IN' not in gdf[g]['jSs'].values and 'N' not in gdf[g]['cigars'].iloc[0]]
            for g in xg:
                del gdf[g]
                
        if len(gdf) >= 1:
            gene = '-'.join(sorted(gdf))
            labels = ['']*len(gdf)
            for i, g in enumerate(sorted(gdf)):
                if len(set(gdf[g]['labels'])) == 1:
                    labels[i] = gdf[g]['labels'].iloc[0]
                else: 
                    labels[i] = 'intron'
            label = '-'.join(labels)
            
            # based on the final gene dataframe (gdf), 
            # a gene name has now been assigned
            
        else:
            
            gene = ''; label = ''
            
    return label, gene, df


###############################################################################
# main
    
# handle mates if they're given
if 

# count reads and assign
cnt = {}
countLabels = set()
#otc = open(output + '-check_assignments_singleMappers.txt', 'w')
f_decisions = open(output+"_singlemapper_decisions.tsv", "w")
with open(bedsingle) as f:
    for i, line in enumerate(f):    
        
        ch, x0, x1, name, strand, gene, info, tlen, cov = line.rstrip().rsplit('\t')
        
        if 'tRNA' in gene:
            
            gene = gene.replace('-','.')
            gene = gene + '_tRNA'
            
        if i == 0:
            
            genes = ["_".join(gene.rsplit("_")[:-1])]; labels = [gene.rsplit("_")[-1]]; infos = [info]; covs = [cov]; tlens = [tlen]
            r0 = name
            
        else:
            
            if name != r0: 
                
                cell, umi = get_cell_UMI(r0, protocol)
                label, xgene, df = gene_assignment_single(genes, labels, infos, covs, tlens)
                
#                if len(df)>1: 
#                    otc.write(r0)
#                    otc.write(str(df))
#                    otc.write('\n' + label + ' '+ xgene)
#                    otc.write('\n\n')
                
                if xgene != '': 
                    cnt = addReadCount(cnt, cell, xgene, umi, label)
                    countLabels.add(label)
                    
                # write mapping to file, for later use in targeted sequencing analysis
                # (bit redundant for single mappers, but easier for later workflow)                    
                f_decisions.write(r0+'\t'+xgene+'\n')

                    
                # start again
                genes = ["_".join(gene.rsplit("_")[:-1])]; labels = [gene.rsplit("_")[-1]]; infos = [info]; covs = [cov]; tlens = [tlen]
                r0 = name
                
            else: 
                
                genes.append("_".join(gene.rsplit("_")[:-1])); labels.append(gene.rsplit("_")[-1]); infos.append(info); covs.append(cov); tlens.append(tlen)


#otc.close()
f_decisions.close()

# again with multi-mappers
#otc = open(output + '-check_assignments_multipleMappers.txt', 'w')
f_decisions = open(output+"_multimapper_decisions.tsv", "w")
with open(bedmulti) as f:
    for i, line in enumerate(f):        
        
        # the bed files will contain multiple lines for each of the mappings.
        # to which read they belong can be identified by the read name
        # so per name we collect the (multiple) mappings belonging to that name
        
        ch, x0, x1, name, strand, gene, info, tlen, cov = line.rstrip().rsplit('\t')
        
        if 'tRNA' in gene:
            
            gene = gene.replace('-','.')
            gene = gene + '_tRNA'
            
        if i == 0:
            
            genes = ["_".join(gene.rsplit("_")[:-1])]; labels = [gene.rsplit("_")[-1]]; infos = [info]; covs = [cov]; tlens = [tlen]
            r0 = name # r0 contains the name we're currently working on
            
        else:
            
            # if entry relates to other name (ie read)
            # handle currently collected entries
            if name != r0: 
                                
                cell, umi = get_cell_UMI(r0, protocol)
                label, xgene, df = gene_assignment_single(genes, labels, infos, covs, tlens)
#                if len(df) > 1:
#                    otc.write(r0)
#                    otc.write(str(df))
#                    otc.write('\n' + label + ' '+ xgene)
#                    otc.write('\n\n')
                
                # write decision to file, for later use in targeted sequencing analysis
                f_decisions.write(r0+'\t'+xgene+'\n')
                
                if xgene != '':
                    cnt = addReadCount(cnt, cell, xgene, umi, label)
                    countLabels.add(label)
                    
                # start again
                genes = ["_".join(gene.rsplit("_")[:-1])]; labels = [gene.rsplit("_")[-1]]; infos = [info]; covs = [cov]; tlens = [tlen]
                r0 = name
                
            else: 
                
                # otherwise just keep collecting entries 
                genes.append("_".join(gene.rsplit("_")[:-1])); labels.append(gene.rsplit("_")[-1]); infos.append(info); covs.append(cov); tlens.append(tlen)
#otc.close()
f_decisions.close()

###############################################################################
# more functions
                
cntdf = pd.DataFrame(cnt)
del cnt

pickle.dump(cntdf, open(output + '.pickle', 'wb'))

#
def sumCounts(x1, x2):
    if type(x1) == dict and type(x2) == dict:
        x = x1
        for umi in x2:
            if umi not in x1:
                x[umi] = x2[umi]
            else:
                x[umi] = x[umi]+x2[umi]
    elif type(x1) == dict and type(x2) != dict:
        x = x1
    elif type(x1) != dict and type(x2) == dict:
        x = x2
    else:
        x = x1
    return x

def aggCounts(xs):
    ax = xs[0]
    if any([type(x) == dict for x in xs]):
        for x in xs[1:]:
            ax = sumCounts(ax,x)
    return ax

# count tables
# total reads
def countTotalReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']: 
        umi = 'A'
        y = sum(x[umi].values()) if type(x) == dict else 0
    return y

def countTotalUMI(x):
    return len(x) if type(x) == dict else 0

# reads with no introns
def countExonReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x if 'intron' not in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']:
        umi = 'A'
        y = 0
        if type(x) == dict:
            y = sum([x[umi][k] if 'exon' in k else 0 for k in x[umi]])
    return y

def countExonUMI(x):
    return len([u for u in x if 'intron' not in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0

# unspliced reads
def countIntronReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x', 'smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x if 'intron' in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']: 
        umi = 'A'
        y = 0
        if type(x) == dict:
            y = sum([x[umi][k] if 'intron' in k else 0 for k in x[umi]])
    return y

def countIntronUMI(x):
    return len([u for u in x if 'intron' in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0

# estimated transcripts 
K = 4**len(umi)
def bc2trans(x):
    if x >= K:
        t = np.log(1.-(float(K)-1e-3)/K)/np.log(1.-1./K)
    elif x > 0 and x < K:
        t = np.log(1.-float(x)/K)/np.log(1.-1./K)
    elif x == 0:
        t = 0
    return  t

###############################################################################
# more main script
    
fout = open(output + '_mapStats.log', 'w')

total_reads_df = cntdf.applymap(lambda x: countTotalReads(x, protocol))
total_umi_df = cntdf.applymap(lambda x: countTotalUMI(x))
transcripts_total_df = total_umi_df.applymap(bc2trans)

fout.write('Total mapped reads:\t' + str(total_reads_df.sum().sum())  + '\n')
fout.write('Dimension of raw dataset:\t' + str(cntdf.shape) + '\n')

total_reads_df.to_csv(output + '_total.ReadCounts.tsv', sep = '\t')
total_umi_df.to_csv(output + '_total.UFICounts.tsv', sep = '\t')
transcripts_total_df.to_csv(output + '_total.TranscriptCounts.tsv', sep = '\t')

# select tRNA and genes indexes:
tRNAs = [idx for idx in cntdf.index if 'tRNA' in idx]
genes = [idx for idx in cntdf.index if idx not in tRNAs]

fout.write('Total reads assigned to tRNA:\t' + str(total_reads_df.loc[tRNAs].sum().sum())  + '\n')
fout.write('Number of uni/multi-tRNA detected:\t' + str(len(tRNAs))  + '\n')

# tRNA tables
cntdf_tRNA = cntdf.loc[tRNAs]
total_reads_tRNA = cntdf_tRNA.applymap(lambda x: countTotalReads(x, protocol))
total_UMI_tRNA = cntdf_tRNA.applymap(countTotalUMI)

total_reads_tRNA['type'] = ['-'.join(sorted(set([t.rsplit('.')[-1] for t in idx.rsplit('-')]))) for idx in total_reads_tRNA.index]
total_reads_tRNA = total_reads_tRNA.groupby('type').aggregate(sum)
total_UMI_tRNA['type'] = ['-'.join(sorted(set([t.rsplit('.')[-1] for t in idx.rsplit('-')]))) for idx in total_UMI_tRNA.index]
total_UMI_tRNA = total_UMI_tRNA.groupby('type').aggregate(sum)

total_reads_tRNA.to_csv(output + '_tRNA.ReadCounts.tsv', sep = '\t')
total_UMI_tRNA.to_csv(output + '_tRNA.UFICounts.tsv', sep = '\t')

fout.write('Number of tRNA detected after collapsing:\t' + str(len(total_reads_tRNA))  + '\n')

# gene tables
uni_genes = [g for g in genes if '-' not in g]
fout.write('Total reads assigned to genes:\t' + str(total_reads_df.loc[genes].sum().sum())  + '\n')
fout.write('Number of uni/multi-genes detected:\t' + str(len(genes))  + '\n')
fout.write('Total reads assigned to uni-genes:\t' + str(total_reads_df.loc[uni_genes].sum().sum())  + '\n')
fout.write('Number of uni-genes:\t' + str(len(uni_genes))  + '\n')

cntdf_genes = cntdf.loc[genes].copy()

###############################################################################
# more functions

def reduceGeneName(gene, uni_genes = uni_genes):
    rg = gene
    if gene.count('-') == 0:
        rg = gene
    else:
        if sum([g in uni_genes for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g in uni_genes][0]
        elif sum([g.rsplit('_')[1][:2]!="Gm" for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g.rsplit('_')[1][:2]!="Gm"][0]
    return rg

def fixGeneLabels(xdf):
    for idx in xdf.index:
        if idx != xdf.loc[idx,'new_gene']:
            i = np.argmax([x==xdf.loc[idx,'new_gene'] for x in idx.rsplit('-')])
            for cell in xdf.columns[:-1]:
                if type(xdf.loc[idx,cell]) == dict:
                    for umi in xdf.loc[idx,cell]:
                        xdf.loc[idx,cell][umi] = Counter([k.rsplit('-')[i] for k in xdf.loc[idx,cell][umi].elements()])
    return xdf

###############################################################################
# more main

cntdf_genes['new_gene'] = [reduceGeneName(idx) for idx in cntdf_genes.index]
cntdf_genes = fixGeneLabels(cntdf_genes)
agg_cntdf_genes = cntdf_genes.groupby('new_gene').aggregate(aggCounts)

uni_genes = [g for g in agg_cntdf_genes.index if '-' not in g]
multi_genes = [g for g in agg_cntdf_genes.index if '-' in g]

total_aggreads_df = agg_cntdf_genes.applymap(lambda x: countTotalReads(x))

fout.write('Total reads after gene aggregation:\t' + str(total_aggreads_df.sum().sum())  + '\n')
fout.write('Dimension of gene-dataset after aggregation:\t' + str(total_aggreads_df.shape) + '\n')
fout.write('Total reads assigned to uni-gene after aggregation:\t' + str(total_aggreads_df.loc[uni_genes].sum().sum())  + '\n')
fout.write('Number of uni-genes after aggregation:\t' + str(len(uni_genes)) + '\n')
fout.write('Total reads assigned to multi-gene after aggregation:\t' + str(total_aggreads_df.loc[multi_genes].sum().sum())  + '\n')
fout.write('Number of multi-gene after aggregation:\t' + str(len(multi_genes))  + '\n')

multi_genes_singleLabel = []
multi_genes_multiLabel = []
multi_genes_multiTags = []
for g in multi_genes:
    ks = set(['&'.join(sorted(x[umi].keys())) for x in agg_cntdf_genes.loc[g] if type(x) == dict for umi in x])
    if len(ks) == 1:
        ks = set(list(ks)[0].rsplit('-'))
        if len(ks) == 1:
            multi_genes_singleLabel.append(g)
        else:
            multi_genes_multiLabel.append(g)
    else:
        multi_genes_multiTags.append(g)

fout.write('Total reads assigned to multi-genes after aggregation that have single labels (all exons/introns):\t' + str(total_aggreads_df.loc[multi_genes_singleLabel].sum().sum())  + '\n')
fout.write('Number of multi-gene after aggregation:\t' + str(len(multi_genes_singleLabel))  + '\n')
fout.write('Total reads assigned to multi-genes after aggregation that have multiple labels (exon-intron etc):\t' + str(total_aggreads_df.loc[multi_genes_multiLabel].sum().sum())  + '\n')
fout.write('Number of multi-gene after aggregation that have multiple labels:\t' + str(len(multi_genes_multiLabel))  + '\n')
fout.write('Total reads assigned to multi-genes after aggregation that have distinct tags (exon-intron and intron-intron):\t' + str(total_aggreads_df.loc[multi_genes_multiTags].sum().sum())  + '\n')
fout.write('Number of multi-genes after aggregation that have distinct tags (exon-intron and intron-intron):\t' + str(len(multi_genes_multiTags))  + '\n')

fout.close()
uni_genes += multi_genes_singleLabel

# unique gene tables
cntdf_unigenes = agg_cntdf_genes.loc[uni_genes]
cntdf_multigenes = agg_cntdf_genes.loc[[idx for idx in agg_cntdf_genes.index if idx not in uni_genes]]

total_reads_multigenes = cntdf_multigenes.applymap(lambda x: countTotalReads(x, protocol))
total_UMI_multigenes = cntdf_multigenes.applymap(countTotalUMI)
total_transcripts_multigenes = total_UMI_multigenes.applymap(bc2trans)

total_reads_unigenes = cntdf_unigenes.applymap(lambda x: countTotalReads(x, protocol))
spliced_reads_unigenes = cntdf_unigenes.applymap(lambda x: countExonReads(x, protocol))
unspliced_reads_unigenes = cntdf_unigenes.applymap(lambda x: countIntronReads(x, protocol))
total_UMI_unigenes = cntdf_unigenes.applymap(countTotalUMI)
spliced_UMI_unigenes = cntdf_unigenes.applymap(countExonUMI)
unspliced_UMI_unigenes = cntdf_unigenes.applymap(countIntronUMI)
total_transcripts_unigenes = total_UMI_unigenes.applymap(bc2trans)
spliced_transcripts_unigenes = spliced_UMI_unigenes.applymap(bc2trans)
unspliced_transcripts_unigenes = unspliced_UMI_unigenes.applymap(bc2trans)

total_reads_unigenes.to_csv(output + '_uniaggGenes_total.ReadCounts.tsv', sep = '\t')
total_UMI_unigenes.to_csv(output + '_uniaggGenes_total.UFICounts.tsv', sep = '\t')
total_transcripts_unigenes.to_csv(output + '_uniaggGenes_total.TranscriptCounts.tsv', sep = '\t')

total_reads_multigenes.to_csv(output + '_multiaggGenes_total.ReadCounts.tsv', sep = '\t')
total_UMI_multigenes.to_csv(output + '_multiaggGenes_total.UFICounts.tsv', sep = '\t')
total_transcripts_multigenes.to_csv(output + '_multiaggGenes_total.TranscriptCounts.tsv', sep = '\t')

unspliced_reads_unigenes.to_csv(output + '_uniaggGenes_unspliced.ReadCounts.tsv', sep = '\t')
unspliced_UMI_unigenes.to_csv(output + '_uniaggGenes_unspliced.UFICounts.tsv', sep = '\t')
unspliced_transcripts_unigenes.to_csv(output + '_uniaggGenes_unspliced.TranscriptCounts.tsv', sep = '\t')

spliced_reads_unigenes.to_csv(output + '_uniaggGenes_spliced.ReadCounts.tsv', sep = '\t')
spliced_UMI_unigenes.to_csv(output + '_uniaggGenes_spliced.UFICounts.tsv', sep = '\t')
spliced_transcripts_unigenes.to_csv(output + '_uniaggGenes_spliced.TranscriptCounts.tsv', sep = '\t')




