#!/bin/bash


# Now perform STAR mapping and counting of reads.
# Note that featureCounts needs to be installed (but doesn't require path).

# Note that STAR requires quite some memory, e.g. 50G

# Paths
path2star=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin

# Parameters to set (see also below):
"
nrthreads=1
genome=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/STAR-indexed-49
myfastqfile=head_clone-12-no-dox-rep1_AAC3KKWM5_S19_L001
"

################################################################################
# Mapping with STAR

${path2star}/STAR --runThreadN ${nrthreads} --genomeDir ${genome} \
                  --readFilesIn ${myfastqfile}_R1_val_1.fq.gz ${myfastqfile}_R2_val_2.fq.gz \
                  --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within \
                  --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix ${myfastqfile}_

# Note: to view this output, use:
# samtools view -h ${myfastqfile}_Aligned.out.bam

################################################################################
# applying featureCounts
# See http://subread.sourceforge.net/featureCounts.html for further information

# Parameters to set:
"
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
"

featureCounts -Q 10 \
              -T ${nrthreads} \
              -a ${gtffile} \
              -o ${myfastqfile}_gene_assignments.tmp \
              -R BAM ${myfastqfile}_Aligned.out.bam 
              # -M
              
              # Note that the "-M" flag enables counting of multi mappers.
              # Default behavior is to ignore multi-mappers.
              # See notes customized_scripts/countTables_umiTools.sh file
              
# Convert the file to a more readable format              
cut -f 1,7 ${myfastqfile}_gene_assignments.tmp | sed 1d > ${myfastqfile}_matrix.txt
 
              
              
              
