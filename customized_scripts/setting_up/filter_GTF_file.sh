
# this code needs to be executed manually
# convenient to use interactive terminal for this at hpc: 
# srun --nodes=1 -c 1 --mem=10G --time=5:00:00 --pty bash -i

gtffile=Homo_sapiens.GRCh38.93
gtffile=Mus_musculus.GRCm39.107

head -1000 ${gtffile}.gtf > head_${gtffile}.gtf

# Trying it on head only
grep -E "(^\\#\\!)|(gene_biotype \"(antisense|protein_coding|lincRNA)\")" head_${gtffile}.gtf > filtered_head_${gtffile}.gtf

# Actual command
grep -E "(^\\#\\!)|(gene_biotype \"(antisense|protein_coding|lincRNA)\")" ${gtffile}.gtf > filteredcustom_${gtffile}.gtf

# sanity check
# new file doesn't have other biotypes
grep -v -E "(gene_biotype \"(antisense|protein_coding|lincRNA)\")" filteredcustom_${gtffile}.gtf
# whereas in the original file, there are still other biotypes
grep -v -E "(gene_biotype \"(antisense|protein_coding|lincRNA)\")" ${gtffile}.gtf | head