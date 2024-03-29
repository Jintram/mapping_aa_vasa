
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

### For the custom written read counting
- **refBED** This is the annotation of the genome, encoded in a customized bed file. This bed file is generated based on a gtf file, see the file: `./setting_up/convert_gtf_to_bed.sh`.

### For the UMItools read counting
- **gtffile** Standard gtf file.

### Software
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

### Part A, pre-processing and mapping with STAR
The main bash script can be found is: `L_TS_submit_vasaplate_map_LOCAL.sh`, 
this consists of multiple steps:
1. Extract cell barcodes
2. Trimming
3. Ribosomal mapping and discarding rRNA reads
4. Mapping with STAR

### Part B, creating count tables, ↓option 1↓
5. Custom script that deals with multi-mappers
6. Creating count tables with custom Python scripts
7. (Targeted-sequencing specific stuff, not applicable.)

### Part B, creating count tables, ↓option 2↓
5. Use "umiCounts" to convert the mapped reads to UMI counts. This requires the script 
`countTables_umiTools.sh`. 
This script requires as input a `.bam` file and a `.gtf` file.
Which itself conists of 3 steps:

countTables_umiTools script steps:
1. featureCounts is called to count features.
2. The file is re-arranged to accomodate UMI tools input, bam files are sorted, umi_tools is called to count reads. Optionally UMIs are trimmed (this is sort of a hack to deal with UMIs of uneven lengths, which occured in data I handled previously).
3. Temp files are removed.

# Bulk RNA sequencing

## Mapping CEL-Seq2 generated bulk RNA sequencing

For that purpose, the standard pipeline can be used. 

## Mapping paired end data (without UMIs or barcodes)

For that purpose I wrote some scripts that are now in the bulk sequencing folder
(`./bulk_seq_mapping`). See also https://github.com/Jintram/mapping_aa_vasa/tree/master/bulk_seq_mapping.

# How to install a new genome

## Required files

You need three files:
1. The STAR reference files.
2. The gtf file (gene annotation), also a filtered version if you want only transcriptome or other biotype filtering. Additionally, this needs to be converted to a custom format if you want to use the custom annotation scripts that were written for Vasa, but not if you're using "option B" for creating the count tables, ie umiCounts. The umiCounts method directly uses the .gtf file.
3. A fasta file with all ribosomal RNA sequences (this is only used for recognizing and discarding those reads). 

## Downloading reference files.

Reference genomes can be downloaded from http://ensembl.org. 

### ERCC RNA Spike-In Control Mixes

In case Spike-ins were added, you also need the sequences and annotation for these.
You'll need the ERCC92.fa and ERCC92.gtf files, which can be downloaded following the notes below.

#### Download manual

Through our order number I ended up on Thermofisher website, 
https://www.thermofisher.com/order/catalog/product/4456740
Then looked for manuals and protocols
https://www.thermofisher.com/search/results?query=4456740&persona=DocSupport&type=Product+FAQs&filter=document.result_type_s%3AManuals%20%26%20Protocols&refinementAction=true

Manual:
https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_086340.pdf

#### Download files

This manual refers to: www.appliedbiosystems.com.

If you search that website for ERCC92.fa, you find an ERCC92 zip file:
https://www.thermofisher.com/search/results?query=ERCC92.fa&focusarea=Search%20All
--> "ERCC92.fa & ERCC92.gtf sequence and annotation files (.zip)"
--> https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip

This gives you the ERCC92.fa and ERCC92.gtf file that you need.

### Primary assembly fa reference file
At the time of writing, a convenient overview can be found at:
http://www.ensembl.org/info/data/ftp/index.html/
E.g. for mice, follow the link to:
http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/dna/
and then download the file:
`Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`

### Gene annotation file
Above website refers you to:
http://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/
and then you'll need:
`Mus_musculus.GRCm39.107.gtf.gz`

### Getting the files

Both these files can be downloaded directly at a server by using the commands:

```
wget http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/Mus_musculus.GRCm39.107.gtf.gz
```

See the folder `./setting_up/` for conversion scripts that take the .gtf file as input and output the customized bed format, you'll need the script `convert_gtf_to_bed.sh`.

### The ribosomal RNA file

This is a very straightforward file in fasta format, e.g. for human:
```
>12S_RNR1
AATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGTTCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGC(...)
>16S_RNR2
GCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAACCAAGC(...)
(...)
```

Ribosomal genes are found in repeating patterns and are therefor a bit curious to annotate. E.g. for mice, see also https://www.ncbi.nlm.nih.gov/gene?cmd=retrieve&list_uids=100861531.

The (I think) custom-made file contains the following sequences for mice:
- rRNA-Rn45s; this is the sequence `NR_046233.2` found at https://www.ncbi.nlm.nih.gov/nuccore/NR_046233.2?report=fasta
- rRNA-Rn5s; this is the sequence `NR_030686.1` found at https://www.ncbi.nlm.nih.gov/nuccore/NR_030686.1?report=fasta
- rRNA-12s_16s; this seems to be sequence `V00665.1` at https://www.ncbi.nlm.nih.gov/nuccore/V00665.1?report=fasta
- rRNA-Rn47s; this seems to be sequence `AH002076.2` at https://www.ncbi.nlm.nih.gov/nuccore/AH002076.2?report=fasta

See also VanInsberghe et al. Nature. 2021 http://doi.org/10.1038/s41586-021-03887-4.

### "Installing" the reference for STAR

You'll need to first perform "indexing" before STAR can use the genome. Since this task
also takes computing power, the command is incorporated in the following script:

```./STAR/generate_index_STAR.sh```

Examples how to start this scripts are given in the `./STAR/submit_index_star.sh` script.

You'll need to navigate to the directory with the genomes (e.g. /hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107/ in my case), copy the `generate_index_STAR.sh` file into that directory, and start the job appropriately, e.g.:
```
fafile=./ensembl/Mus_musculus.GRCm39.dna.primary_assembly.fa
gtffile=./ensembl/Mus_musculus.GRCm39.107.gtf.gz
overhang=150

sbatch --export=ALL,fafile="${fafile}",gtffile="${gtffile}",overhang="$overhang" generate_index_STAR.sh 
```
The overhang parameter is simply the max read length minus one (see documentation).
Note that the script will assume the files are zipped, but requires the filenames without .gz extension.

### Filtering the gtf file.

To filter the gtf file, you can use `./setting_up/filter_GTF_file.sh`.

### Combining ERCC and genome files

The .fa for respectively the ERCC and genome can simply be merged, as can the
.gtf files.

Example commands:

```
cat /hpc/hub_oudenaarden/mwehrens/ref/ERCC/ERCC92.fa /hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107/ensembl/Mus_musculus.GRCm39.dna.primary_assembly.fa > /hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107.ERCC/GRCm39.dna.primary_plus_ERCC.fa

cat /hpc/hub_oudenaarden/mwehrens/ref/ERCC/ERCC92.gtf /hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107/ensembl/filteredcustom_Mus_musculus.GRCm39.107.gtf > /hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107.ERCC/GRCm39.dna.primary_plus_ERCC.gtf
```







