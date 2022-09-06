
# Mapping pipeline for RNA-sequencing data

## Acknowledgements
These scripts are based on Anna Alemany's Vasa-seq mapping pipeline.
This pipeline can be found at https://github.com/hemberg-lab/VASAseq_2022.
See also: 
- Salmen F, Jonghe J De, Kaminski TS, Alemany A, Parada GE, Verity-legg J, Yanagida A, Kohler TN, Battich N, Brekel F Van Den, Ellermann AL, Arias AM, Nichols J, Hemberg M, Hollfelder F, Oudenaarden A Van. High-throughput total RNA sequencing in single cells using VASA-seq. Nat Biotechnol. Published online 2022. http://doi.org/10.1038/s41587-022-01361-8.

Scripts customized and adapted by Martijn Wehrens. 
The scripts were also used for (among others) the following publication(s):
- Wehrens M, Leeuw AE de, Wright-Clark M, Eding JEC, Boogerd CJ, Molenaar B, Kraak PH van der, Kuster DWD, Velden J van der, Michels M, Vink A, Rooij E van. Single-cell transcriptomics provides insights into hypertrophic cardiomyopathy. Cell Rep. 2022;39(6). http://doi.org/10.1016/j.celrep.2022.110809

## Motivation
I needed a mapping pipeline to deal with RNA-sequencing data from diverse
sources, and did not want to start from scratch. Anna Alemany / Single Cell Discoveries
where kind enough to share their code, now available at repo listed above.

This mapping pipeline should also be able to (as much as possible) reproduce what the Cell Ranger 
scripts from 10X are doing. (See also methods Wehrens 2022).

This is an ongoing project, and many things in this repo still need to be revised (eg filenames that don't make sense, etc).

## Requirements
(To add; STAR etc.)

## Setting up
References need to set up in a custom way.
Reference files that are needed are:
- **riboref** This is a fasta format file that contains all the rRNA sequences.
- **genome** This is the STAR indexed genome.
- **refBED** This is the annotation of the genome, encoded in a customized bed file. This bed file is generated based on a gtf file, see the file: `/mapping_aa_private_vasa/customized_scripts/setting_up/convert_gtf_to_bed.sh`.

## File structure and usage
..







