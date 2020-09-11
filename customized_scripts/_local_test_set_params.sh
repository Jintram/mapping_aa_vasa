# For running local
p2s=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts
p2trimgalore=/Users/m.wehrens/Software_custom/TrimGalore-0.6.6
path2trimgalore=/Users/m.wehrens/Software_custom/TrimGalore-0.6.6
p2cutadapt=/Users/m.wehrens/Library/Python/3.7/bin
path2cutadapt=/Users/m.wehrens/Library/Python/3.7/bin
lib=MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_L001
file2trim=${lib}_TS
file2trim=${lib}_pT_R2
file2trim=${lib}_nc_R2
file2trim=${lib}
cutoff=10

p2bwa=/Users/m.wehrens/Software_custom/bwa
fq_R1=${lib}_TS_R1_cbc_val_1_HATCG.fq.gz
fq_R2=${lib}_TS_R2_cbc_val_2_HATCG.fq.gz
riboref=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/unique_rRNA_human.fa
ref=$riboref
#p2samtools=/Users/m.wehrens/Software_custom/samtools-1.9_bin/bin
p2samtools=/Users/m.wehrens/Software_custom/samtools-1.2
out=$lib