# looking at data sjoerd

zcat myfile.zip | head -n 10 

# Sanity check, checking structure of files
# indeed "IX" files just contain the index sequence
f=clone-11-dox-rep1_AAC3KKWM5_S10_L001_I1_001.fastq.gz
zcat $f | head -n 30 

# These should hold sequences:
f=cone-5-dox-rep1_AAC3KKWM5_S22_L001_R1_001.fastq.gz
zcat $f | head -n 30 
  # note: sequences are 51 nt long
  
