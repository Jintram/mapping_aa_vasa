# looking at data sjoerd
# 2022-09-06

# This file describes what I did and some additional stuff to look at the files. 

# "Simple" bulk RNA seq data, already demultiplexed, no UMIs involved, index
# primers already removed.

zcat myfile.zip | head -n 10 

# Sanity check, checking structure of files
# indeed "IX" files just contain the index sequence
f=clone-11-dox-rep1_AAC3KKWM5_S10_L001_I1_001.fastq.gz
zcat $f | head -n 30 

# These should hold sequences:
f=cone-5-dox-rep1_AAC3KKWM5_S22_L001_R1_001.fastq.gz
zcat $f | head -n 30 
  # note: R1 sequences are 51 nt long
f=cone-5-dox-rep1_AAC3KKWM5_S22_L001_R2_001.fastq.gz
zcat $f | head -n 30 
  
# Note that all files come from lane 1 (L001).

# File structure is like:
# c(l)one-XXX-ZZZ-repYYY
# with XXX the clone # and YYY the rep, ZZZ gives "dox" or "no-dox" condition
# there's one typo "cone" instead of "clone"
# The samples are:
# ls -lha *_R2_* | awk '{print $9}' | sed 's/_L001_R2_001.fastq.gz//g'
clone-11-dox-rep1_AAC3KKWM5_S10
clone-11-dox-rep2_AAC3KKWM5_S11
clone-11-dox-rep3_AAC3KKWM5_S12
clone-11-no-dox-rep1_AAC3KKWM5_S13
clone-11-no-dox-rep2_AAC3KKWM5_S14
clone-11-no-dox-rep3_AAC3KKWM5_S15
clone-12-dox-rep1_AAC3KKWM5_S16
clone-12-dox-rep2_AAC3KKWM5_S17
clone-12-dox-rep3_AAC3KKWM5_S18
clone-12-no-dox-rep1_AAC3KKWM5_S19
clone-12-no-dox-rep2_AAC3KKWM5_S20
clone-12-no-dox-rep3_AAC3KKWM5_S21
clone-2-dox-rep1_AAC3KKWM5_S1
clone-2-dox-rep2_AAC3KKWM5_S2
clone-2-dox-rep3_AAC3KKWM5_S3
clone-3-dox-rep1_AAC3KKWM5_S4
clone-3-dox-rep2_AAC3KKWM5_S5
clone-3-dox-rep3_AAC3KKWM5_S6
clone-5-no-dox-rep1_AAC3KKWM5_S7
clone-5-no-dox-rep2_AAC3KKWM5_S8
clone-5-no-dox-rep3_AAC3KKWM5_S9
cone-5-dox-rep1_AAC3KKWM5_S22
cone-5-dox-rep2_AAC3KKWM5_S23
cone-5-dox-rep3_AAC3KKWM5_S24
# (24 samples)

# Let's create a test dataset.
f1=clone-12-no-dox-rep1_AAC3KKWM5_S19_L001_R1_001.fastq
zcat $f1.gz | head -n 1000 > head_$f1
f2=clone-12-no-dox-rep1_AAC3KKWM5_S19_L001_R2_001.fastq
zcat $f2.gz | head -n 1000 > head_$f2
gzip head_$f1
gzip head_$f2

# Need to remove the suffixes from the filenames
for file in *.gz
do 
    echo "${file%_001.fastq.gz}"
    mv "$file" "${file%_001.fastq.gz}".fast.gz
done

for file in *.gz
do 
    mv "$file" "${file%.fast.gz}".fastq.gz
    #echo "${file%.fast.gz}".fastq.gz
done

# Let's see if we can map the test dataset using my pipline
# Customized the scripts:
# ./projects/Sjoerd_202209/run_parameters_Sjoerd.sh
# Using:
# ./projects/Sjoerd_202209/submit_jobs_test_Sjoerd.sh
# ./projects/Sjoerd_202209/submit_Counttables_umiCounts_Sjoerd_discardmultiandfiltered.sh
#
# --> Doesn't work; of course, that pipeline doesn't handle paired-end
# data. 

# The mapping however is actually quite simple; so let's create a small
# bulk seq pipeline.






















