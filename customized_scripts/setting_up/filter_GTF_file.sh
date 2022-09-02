
head Homo_sapiens.GRCh38.93.gtf > head_Homo_sapiens.GRCh38.93.gtf

# Trying it on head only
grep -E "(^\\#\\!)|(gene_biotype \"(antisense|protein_coding|lincRNA)\")" head_Homo_sapiens.GRCh38.93.gtf > filtered_head_Homo_sapiens.GRCh38.93.gtf

# Actual command
grep -E "(^\\#\\!)|(gene_biotype \"(antisense|protein_coding|lincRNA)\")" Homo_sapiens.GRCh38.93.gtf > filteredcustom_Homo_sapiens.GRCh38.93.gtf

# sanity check
grep -v -E "(gene_biotype \"(antisense|protein_coding|lincRNA)\")" filteredcustom_Homo_sapiens.GRCh38.93.gtf
