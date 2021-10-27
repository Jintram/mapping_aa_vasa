

# umi length count
awk 'BEGIN{FS="\t"} {
  if (substr($0,1,1)=="@") {} 
  else {
    split($1,first_part,";"); 
    print first_part[5] "\t" length(first_part[5])
  }}' ${bamfile}.bam.featureCounts.sam > umi_lengths.txt