#!/bin/bash


# Now perform STAR mapping and counting of reads.
# Note that featureCounts needs to be installed (but doesn't require path).

# Note that STAR requires quite some memory, e.g. 50G

# Paths
path2star=/hpc/hub_oudenaarden/mwehrens/bin/miniconda3/bin

# Parameters to set (see also below):
#nrthreads=1
#genome=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/STAR-indexed-49
#myfastqfile=head_clone-12-no-dox-rep1_AAC3KKWM5_S19_L001
#parentdir=/hpc/hub_oudenaarden/mwehrens/fastq/Sjoerd


################################################################################
# Mapping with STAR

cd $TMPDIR

${path2star}/STAR --runThreadN ${nrthreads} --genomeDir ${genome} \
                  --readFilesIn ${parentdir}/workdir/${myfastqfile}_R1_val_1.fq.gz ${parentdir}/workdir/${myfastqfile}_R2_val_2.fq.gz \
                  --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within \
                  --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix ${myfastqfile}_

ls $TMPDIR/${myfastqfile}_Aligned.out.bam

mv $TMPDIR/${myfastqfile}_Aligned.out.bam ${parentdir}/workdir/
cd ${parentdir}/workdir/

# Note: to view this output, use:
# samtools view -h ${myfastqfile}_Aligned.out.bam

################################################################################
# applying featureCounts
# See http://subread.sourceforge.net/featureCounts.html for further information

# Parameters to set:
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/filteredcustom_Homo_sapiens.GRCh38.93.gtf
  # only protein_coding, lincRNA, antisense (like cell ranger)

featureCounts -Q 10 \
              -T ${nrthreads} \
              -a ${gtffile} \
              -o ${parentdir}/workdir/${myfastqfile}_gene_assignments.tmp \
              -R BAM ${parentdir}/workdir/${myfastqfile}_Aligned.out.bam
              # -M
              
              # Note that the "-M" flag enables counting of multi mappers.
              # Default behavior is to ignore multi-mappers.
              # See notes customized_scripts/countTables_umiTools.sh file
              
# Convert the file to a more readable format              
cut -f 1,7 ${parentdir}/workdir/${myfastqfile}_gene_assignments.tmp | sed 1d > ${parentdir}/workdir/${myfastqfile}_matrix.txt

mv ${parentdir}/workdir/${myfastqfile}_matrix.txt ${parentdir}/counttables/${myfastqfile}_matrix.txt
 
echo "Part 2 done for ${myfastqfile}.."
              
              
              
