
# Following:
# SubreadUsersGuide.pdf
#
# Not sure what they do with multi-mappers though; it seems
# they also either recommend either not counting them or
# counting them towards both genes..


library(Rsubread)

featureCounts_out2 = 
    featureCounts(  files="/Volumes/fastq_m.wehrens/Mapping/WANG4/mapping/head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.bam",
                    annot.ext="/Volumes/fastq_m.wehrens/Mapping/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf",
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    #GTF.attrType="gene_id",
                    GTF.attrType="gene_name",
                    nthreads = 1)


# So this works rather well;
# Would be easy to write this on a per-cell basis (filter and make little temp files); 
# only question then is how to count UMIs


featureCounts_out2$counts['MTND4P12',]/featureCounts_out2$counts['MT-ND4',]


