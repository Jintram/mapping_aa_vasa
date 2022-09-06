

# New approach, let's try to use RSEM to annotate my reads..
# Uisng https://databeauty.com/blog/tutorial/2016/09/13/RNA-seq-analysis.html


################################################################################
datadir=/Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/testrun3/rsem_test
p2rsem=/Users/m.wehrens/Software_custom/RSEM
p2gtf=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/original_fa_files_human_ref/Homo_sapiens.GRCh38.99.gtf
p2genome=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/original_fa_files_human_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
p2rsemgenome=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/RSEM_prepared/rsem
p2stargenome=/Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/ref_genome_anna/star_v273a_NOMASK_NOERCC_index_74
samplename=MW-TS-S1-Hybr-NB_HTMH2BGXF_S23
readin1=${samplename}_cat_TS_cbc_val_HATCG_R1.nonRibo.fastq
readin2=${samplename}_cat_TS_cbc_val_HATCG_R2.nonRibo.fastq
################################################################################

# prepare a rsem genome first 
cd $p2rsemgenome
${p2rsem}/rsem-prepare-reference --gtf $p2gtf $p2genome rsem
     
# See what happens when we use it on test:
#THIS DOESNT WORK because RSEM requires mapping to transcriptome, but this can be done by STAR easily
#sample=MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo
#${p2rsem}/rsem-calculate-expression --bam --no-bam-output -p 20 --paired-end --forward-prob 1 /Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/testrun3/rsem_test/MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.bam /Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/RSEM_prepared/rsem output/test/rsem 

# First map to transriptome:
STAR --genomeDir $p2stargenome --sjdbGTFfile $p2gtf \
    --readFilesIn $readin1 $readin2 \
    --readFilesCommand cat --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 \
    --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic \
    --runThreadN 20 --outFileNamePrefix star_${sample}_ --sjdbOverhang 74 \
    --outSAMattrRGline ID:$sample  \
    --genomeSAindexNbases 10 --genomeSAsparseD 2 # these two result in lower speed and memory usage (allows running locally)

# Convert to sam and also remove reads that don't have both mates mapped
samtools view -f 2 -h -o star_${sample}_Aligned.toTranscriptome.out.f2.sam star_${sample}_Aligned.toTranscriptome.out.bam

# Now sort out per barcode
$pythonbin ${p2s}/mw_rsem_split_mapped_bam_per_bc.py ${inbam%.bam}.f2.sam ${inbam%.bam}.singlemappers.sam ${inbam%.bam}.multimappers.sam

# GREP for MYBPC3 reads:
# grep "ENST00000545968\|ENST00000256993\|ENST00000399249\|ENST00000544791" star_MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_Aligned.toTranscriptome.out.sam

# TODO?
echo "perhaps add \"-f 2\" filter also to remove one-read mappings?"

# calculate expression  
${p2rsem}/rsem-calculate-expression --seed-length 5 --bam --no-bam-output -p 20 --paired-end --forward-prob 1 \
     /Users/m.wehrens/Data_notbacked/mapping-test-local/pilot3/testrun3/rsem_test/star_MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_Aligned.toTranscriptome.out.bam \
     /Volumes/workdrive_m.wehrens_hubrecht/reference_genomes/RSEM_prepared/rsem output/test/rsem >& \
    output/test/rsem.log
    # Seed length is set to 5 to allow for the shorter R1 reads (with valuable
    # information) to be included.
     
  
# calculate expression per barcode   
while read F  ; do
  
  ${p2rsem}/rsem-calculate-expression --seed-length 5 --bam --no-bam-output -p 20 --paired-end --forward-prob 1 \
       ${datadir}/percell/${F}.sam \
       ${p2rsemgenome} percell/rsem_${F} >& \
      percell/rsem_${F}.log
      
done < ./percell/cells_with_reads.tsv
 
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     