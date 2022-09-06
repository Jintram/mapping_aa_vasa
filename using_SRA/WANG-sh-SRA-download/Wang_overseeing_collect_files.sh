#!/bin/bash

# collect files separately


# The disadvantage of doing all separately is that this will take up
# large amounts of space ...
#for IDENTIFIER in p.N1.plate.97474.part.1 p.N1.plate.97474.part.2 p.N2.plate.97452.part.1 p.N2.plate.97452.part.2 p.N2.plate.97493.part.1 p.N2.plate.97493.part.2 p.N3.plate.97438.part.1 p.N3.plate.97438.part.2 p.N4.plate.97461.part.1 p.N4.plate.97461.part.2 p.N5.plate.97458.part.1 p.N5.plate.97458.part.2 p.N13.plate.100355.part.1 p.N13.plate.100355.part.2 p.N14.plate.104720.part.1 p.N14.plate.104720.part.2
#do
#  sbatch --output=slurm-${sample}-%x.%j.out --job-name=${sample} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${IDENTIFIER}" Wang_collect_files.sh 
#done

# So let's do it in batches of two instead

############################################################
# Also, let's first submit these jobs

samples="p.N1.plate.97474.part.1 p.N1.plate.97474.part.2"
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 

samples="p.N2.plate.97452.part.1 p.N2.plate.97452.part.2"
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 

samples="p.N2.plate.97493.part.1 p.N2.plate.97493.part.2" 
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 

samples="p.N3.plate.97438.part.1 p.N3.plate.97438.part.2" 
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 

############################################################
# And these slightly later

samples="p.N4.plate.97461.part.1 p.N4.plate.97461.part.2" 
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 

samples="p.N5.plate.97458.part.1 p.N5.plate.97458.part.2" 
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 

samples="p.N13.plate.100355.part.1 p.N13.plate.100355.part.2" 
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 

samples="p.N14.plate.104720.part.1 p.N14.plate.104720.part.2"
jname=$(echo "$samples" | sed -r 's/ /_/g')
sbatch -c 2 --output=slurm-${jname}-%x.%j.out --job-name=${jname} --time=48:00:00  --mem=30G --export=ALL,INPUT_IDENTIFIERS="${samples}" Wang_collect_files.sh 
 


