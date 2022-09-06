
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
- **refBED** This is the annotation of the genome, encoded in a customized bed file. This bed file is generated based on a gtf file, see the file: `./setting_up/convert_gtf_to_bed.sh`.

There's a file called `general_parameters.sh` which contains file paths to the
software that needs to be installed.
- path to directory with these custom scripts
- path to trimgalore
- path to cutadapt
- path to bwa
- path to samtools
- path to STAR
- path to python
- path to bedtools

## Input files
This mapping pipeline requires fastq files, archived (.gz format), and it expects 
files from multiple lanes to be merged.

## File structure and usage
Practically, the script you need to edit beforehand to convey the right settings
is 
`run_parameters.sh` that needs to be made specifically for each of your data sets.
In this script, you'll convey the output directory, the read length (minus 1), 
reference files to be used, the protocol, whether targeted-sequencing analysis 
should be performed (no unless you know what this means), and whether the protocol 
is stranded (usually yes).

The main bash script can be found is: `L_TS_submit_vasaplate_map_LOCAL.sh`, 
this consists of multiple steps:
1. Extract cell barcodes
2. Trimming
3. Ribosomal mapping and discarding rRNA reads
4. Mapping with STAR
5. Custom script that deals with multi-mappers
6. Creating count tables with custom Python scripts
7. (Targeted-sequencing specific stuff, not applicable.)







