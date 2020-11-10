

# I think there's something wrong with how the .bed gets filtered into the _genes.bed thing 
# since 

# R1
grep "MYBPC3" MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.f2.singlemappers.R1.intersect.bed | wc -l
# = 129

# R2
grep "MYBPC3" MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.f2.singlemappers.R2.intersect.bed | wc -l
# 431

# So both numbers that are plausible (R2 might cover multiple exons)
grep "MYBPC3" MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.R1.singlemappers_genes.bed | wc -l
# 0 
grep "MYBPC3" MW-TS-S1-Hybr-NB_HTMH2BGXF_S23_cat_TS.nonRibo_E99_Aligned.out.R2.singlemappers_genes.bed | wc -l
# 188

