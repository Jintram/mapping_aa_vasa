#!/bin/bash

# NOTE THAT I'M NOW CURRENTLY LOOKING INTO THE WANG MAPPING, 
# BUT THIS OF COURSE ALSO HAS THESE GENES (OR AT LEAST, IT SHOULD ..)


# First create a test file
mapping_dir=/Volumes/fastq_m.wehrens/Mapping/WANG4/
bigtestfile=p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out
# convert to sam
# samtools view -h -o ${bigtestfile}.sam ${bigtestfile}.bam
# take head
# head -n 100000 ${bigtestfile}.sam > head_${bigtestfile}.sam
# Now merge these two steps
cd $mapping_dir
samtools view -h ${bigtestfile}.bam | head -n 100000 > head_${bigtestfile}.sam

# Convert to bam
testfile=head_${bigtestfile}
samtools view -bS ${testfile}.sam > ${testfile}.bam

# First perhaps do simple manual test, using the python script;
# For this, we first need to generate the overlap bed file
p2s=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts
run_parameter_filepath=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/projects/WANG/run_parameters_Wang_local.sh
general_parameter_filepath=/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/general_parameters_local.sh
lib=head_p.N1.plate.97474.part.1_cat
${p2s}/L_TS_deal_with_singleandmultimappers_paired_stranded.sh $general_parameter_filepath $run_parameter_filepath ${lib}_nc.nonRibo_E99_Aligned.out.bam n

# First, an interesting question is, whether we can find these pseudogenes
# in the single mappers
# Top 1 gene for Rooij cluster 0
grep _MTND4P12_ head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.singlemappers_genes.bed | wc -l
grep _MTND4P12_ head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.multimappers_genes.bed | wc -l  
  # 25 x in the single mappers
  # 239 x in the multi mappers
# Top 2 gene for Rooij cluster 0
grep _MTATP6P1_ head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.singlemappers_genes.bed | wc -l
grep _MTATP6P1_ head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.multimappers_genes.bed | wc -l  
  #  301 in single, 3130 in multi
# Top 3 gene for Rooij cluster 0
grep _MTND1P23_ head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.singlemappers_genes.bed | wc -l
grep _MTND1P23_ head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.multimappers_genes.bed | wc -l  
  # 1x single, 1370x multi

# First, visually inspect what other genes are there with the multis
# open with atom: [head_p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.multimappers_genes.bed] 

# Now, check with Python whether these multi's are called well






