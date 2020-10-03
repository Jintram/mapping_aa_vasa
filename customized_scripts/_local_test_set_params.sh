# For running local
# step1 extract barcodes
p2s=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts
p2trimgalore=/Users/m.wehrens/Software_custom/TrimGalore-0.6.6
path2trimgalore=/Users/m.wehrens/Software_custom/TrimGalore-0.6.6
p2cutadapt=/Users/m.wehrens/Library/Python/3.7/bin
path2cutadapt=/Users/m.wehrens/Library/Python/3.7/bin
lib=MW-TS-S1A-Hybr-Bst_HTMH2BGXF_S25_L001


# TS_extractBC -- needs to be performed manually by command line
pythonbin=/Users/m.wehrens/anaconda3/bin/python
outfq=MW-TS-S1A-Hybr-Bst_HTMH2BGXF_S25_L001
protocol=celseq2
path2scripts=${p2s}
outdir=./

# then TS_trim_paired.sh & TS_trim.sh 
cutoff=10
# For different files:
# paired
file2trim=${lib}_TS
# poly-T
file2trim=${lib}_pT_R2
# non-classified
file2trim=${lib}_nc_R2

# TS_ribo-bwamem.sh and TS_ribo-bwamem_paired.sh
riboref=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/unique_rRNA_human.fa
ref=$riboref
p2bwa=/Users/m.wehrens/Software_custom/bwa
p2samtools=/Users/m.wehrens/Software_custom/samtools-1.2
stranded=y
# for paired
fq_R1=${lib}_TS_R1_cbc_val_1_HATCG.fq.gz
fq_R2=${lib}_TS_R2_cbc_val_2_HATCG.fq.gz
out=${lib}_TS
# for single end (poly-T)
fq=${lib}_pT_R2_cbc_trimmed_HATCG.fq.gz
out=${lib}_pT
# for single end (non-classified)
fq=${lib}_nc_R2_cbc_trimmed_HATCG.fq.gz
out=${lib}_nc

# run TS_map_star_paired.sh and TS_map_star.sh
genome=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/star_v273a_NOMASK_NOERCC_index_74
p2star=/usr/local/bin
# First paired
inputfq1=${lib}_TS_R1.nonRibo.fastq 
inputfq2=${lib}_TS_R2.nonRibo.fastq 
  # locally, it needs to be unzipped first, since zcat doesn't work
outprefix=${lib}_TS.nonRibo.fastq_E99_
# Then the single ones
inputfq=${lib}_pT.nonRibo.fastq
outprefix=${lib}_pT.nonRibo.fastq_E99_

# convert bed files, do some administration for later annotation
refBED=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/Homosapines_ensemble99.homemade_IntronExonTrna.bed
stranded=y
# paired
inbam=${lib}_TS.nonRibo.fastq_E99_Aligned.out.bam
# non-paired
inbam=${lib}_pT.nonRibo.fastq_E99_Aligned.out.bam

# Then create count tables
# use code in TS_submit_..


################################################################################
# manual test on server

lib=MW-TS-S1-Hybr-NB_HTMH2BGXF_S23
ref=HUMAN
n=74
out=XtestX







