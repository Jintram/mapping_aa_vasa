#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 14:17:30 2020

@author: m.wehrens
"""

###############################################################################
# This tests Anna's case where a hit overlaps with an exon-exon junction in
# one gene, and just falls within another gene (at same location)

# minimal test case example
genes, labels, infos, covs, tlens

genes = ['geneA','geneA','geneB']
labels = ['exon','exon','exon']
infos = ['C1:40M;C2:35M5N35M;nM:0;jS:3;XM:2',
         'C1:40M;C2:35M5N35M;nM:0;jS:5;XM:2',
         'C1:40M;C2:35M5N35M;nM:0;jS:IN;XM:2']
covs=[0,0,0]
tlens=[0,0,0]

label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens)
gene_out


# For original version
genes = ['geneA','geneA','geneB']
labels = ['exon','exon','exon']
infos = ['CG:35M50N35M;nM:0;jS:5',
         'CG:35M50N35M;nM:0;jS:3',
         'CG:35M50N35M;nM:0;jS:IN']
covs=[0,0,0]
tlens=[0,0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens)
gene_out

###############################################################################
# Another test case, where a hit overlaps multiple exon-exon junctions of one gene,
# and has only one mismatch with another gene

# For my version
genes = ['geneA','geneA','geneA','geneB']
labels = ['exon','exon','exon','exon']
infos = ['C2:20M10N20M10N20M;nM:0;jS:5;XM:2',
         'C2:20M10N20M10N20M;nM:0;jS:OUT;XM:2',
         'C2:20M10N20M10N20M;nM:0;jS:3;XM:2',
         'C2:30M1N30M;nM:1;jS:IN;XM:2']
covs=[0,0,0,0]
tlens=[0,0,0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens)
gene_out

# For orignial version
genes = ['geneA','geneA','geneA','geneB']
labels = ['exon','exon','exon','exon']
infos = ['CG:20M10N20M10N20M;nM:0;jS:5',
         'CG:20M10N20M10N20M;nM:0;jS:OUT',
         'CG:20M10N20M10N20M;nM:0;jS:3',
         'CG:30M1N30M;nM:1;jS:IN']
covs=[0,0,0,0]
tlens=[0,0,0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens)
gene_out

###############################################################################
# Another test case, hit covers two exons of two separate genes, and has
# appropriate gap 

# For my version
genes = ['geneA','geneB']
labels = ['exon','exon']
infos = ['CG:30M10N30M;nM:0;jS:5;XM:2',
         'CG:30M10N30M;nM:0;jS:3;XM:2']
covs=[0,0]
tlens=[0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens)
gene_out

# Now if this regards mates
genes = ['geneA','geneB']
labels = ['exon','exon']
infos = ['C2:30M;nM:0;jS:IN;XM:2',
         'C1:20M;nM:0;jS:IN;XM:1']
covs=[0,0]
tlens=[0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens)
gene_out
# ---> gives different outcome
# one could pose an additional rule to require both mates be present for a
# gene to remain to be considered


# But.. What now if it hits two exons of geneA? (For my version)
genes = ['geneA','geneA','geneB']
labels = ['exon','exon','exon']
infos = ['CG:30M10N30M10N30M;nM:0;jS:5;XM:2',
         'CG:30M10N30M10N30M;nM:0;jS:3;XM:2',
         'CG:30M10N30M10N30M;nM:0;jS:3;XM:2']
covs=[0,0,0]
tlens=[0,0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens)
gene_out

# And is this compatible with mates? (For my version)
genes = ['geneA','geneA','geneB']
labels = ['exon','exon','exon']
infos = ['C1:20M;nM:0;jS:5;XM:1',
         'C2:30M10N30M10N30M;nM:0;jS:3;XM:2',
         'C2:30M10N30M10N30M;nM:0;jS:3;XM:2']
covs=[0,0,0]
tlens=[0,0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens, bothmates=1)
gene_out

# Test case where double-mate coverage requirement is removed (by "hacking" the call a bit with CG/bothmates)
genes = ['geneA','geneA','geneB']
labels = ['exon','exon','exon']
infos = ['CG:20M;nM:0;jS:5;XM:2',
         'CG:30M10N30M10N30M;nM:0;jS:3;XM:2',
         'CG:30M10N30M10N30M;nM:0;jS:3;XM:2']
covs=[0,0,0]
tlens=[0,0,0]
label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens, bothmates=0)
gene_out


###############################################################################
# "real case": maps to exons one gene and lncRNA of other gene

genes = ['ENSG00000135720_DYNC1LI2_ProteinCoding',
 'ENSG00000135720_DYNC1LI2_ProteinCoding',
 'ENSG00000287965_AC018557.3_lncRNA',
 'ENSG00000135720_DYNC1LI2_ProteinCoding',
 'ENSG00000287965_AC018557.3_lncRNA']
labels = ['exon', 'exon', 'intron', 'exon', 'intron']
infos = ['C1:3M1743N37M;C2:70M1S;nM:0;jS:3;XM:1',
 'C1:3M1743N37M;C2:70M1S;nM:0;jS:5;XM:1',
 'C1:3M1743N37M;C2:70M1S;nM:0;jS:IN;XM:2',
 'C1:3M1743N37M;C2:70M1S;nM:0;jS:IN;XM:2',
 'C1:3M1743N37M;C2:70M1S;nM:0;jS:IN;XM:1']
covs = ['0.434985', '0.434985', '0.416014', '0.565927', '0.586269']
tlens = ['3434', '338', '11727', '338', '11727']

label_out, gene_out, df_out = gene_assignment_single(genes, labels, infos, covs, tlens, bothmates=1)
gene_out











