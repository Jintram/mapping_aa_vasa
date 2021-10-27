#!/bin/sh 

mem=32G
tmpspace=10G

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04/

#gtffile=/Volumes/fastq_m.wehrens/Mapping/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/filteredcustom_Homo_sapiens.GRCh38.93.gtf
  # only protein_coding, lincRNA, antisense (like cell ranger)
  
#nrthreads=
trimUMI=10
#steptorun
#bamfile=..
countmulti=yes

# for testing purposes:
# lib=head_p.N1.plate.97474.part.1_cat

for lib in p.N1.plate.97474.part.1_cat p.N1.plate.97474.part.2_cat p.N2.plate.97452.part.1_cat  p.N2.plate.97452.part.2_cat  p.N2.plate.97493.part.1_cat  p.N2.plate.97493.part.2_cat  p.N3.plate.97438.part.1_cat  p.N3.plate.97438.part.2_cat  p.N4.plate.97461.part.1_cat  p.N4.plate.97461.part.2_cat  p.N5.plate.97458.part.1_cat  p.N5.plate.97458.part.2_cat  p.N13.plate.100355.part.1_cat  p.N13.plate.100355.part.2_cat  p.N14.plate.104720.part.1_cat  p.N14.plate.104720.part.2_cat
do
  echo $lib
  
  bamfile=/hpc/hub_oudenaarden/mwehrens/fastq/WANG4/mapping/${lib}_nc.nonRibo_E99_Aligned.out
  # bamfile=/Volumes/fastq_m.wehrens/Mapping/WANG4/mapping/${lib}_nc.nonRibo_E99_Aligned.out
  
  nrthreads=8
  job1=$(sbatch  --parsable --output=slurm-${lib}-%x.%j.out --job-name=cntUT_${lib} -c ${nrthreads} --gres=tmpspace:${tmpspace} --time=14-00:00:00  --mem=${mem} \
    --export=ALL,gtffile="${gtffile}",nrthreads="${nrthreads}",trimUMI="${trimUMI}",steptorun="1",bamfile="${bamfile}",countmulti="${countmulti}" ${script_path}/countTables_umiTools.sh p.N1.plate.97474.part.1_cat)
  
  nrthreads=1
  job2=$(sbatch --dependency=afterany:${job1} --parsable --output=slurm-${lib}-%x.%j.out --job-name=cntUT_${lib} -c ${nrthreads} --gres=tmpspace:${tmpspace} --time=14-00:00:00  --mem=${mem} \
      --export=ALL,gtffile="${gtffile}",nrthreads="${nrthreads}",trimUMI="${trimUMI}",steptorun="2",bamfile="${bamfile}",countmulti="${countmulti}" ${script_path}/countTables_umiTools.sh p.N1.plate.97474.part.1_cat)

  nrthreads=1
  job3=$(sbatch --dependency=afterany:${job2} --parsable --output=slurm-${lib}-%x.%j.out --job-name=cntUT_${lib} -c ${nrthreads} --gres=tmpspace:${tmpspace} --time=14-00:00:00  --mem=${mem} \
      --export=ALL,gtffile="${gtffile}",nrthreads="${nrthreads}",trimUMI="${trimUMI}",steptorun="3",bamfile="${bamfile}",countmulti="${countmulti}" ${script_path}/countTables_umiTools.sh p.N1.plate.97474.part.1_cat)
    
    
  echo "started job1: $job1; job2: $job2; job3: $job3"  
  
done

# Extra mem was needed for (after initial run with 10G mem):
# for lib in p.N1.plate.97474.part.2_cat p.N14.plate.104720.part.1_cat p.N14.plate.104720.part.2_cat

# For some reason, 32G wasn't enough for p.14, let's just try 128G.
# mem=128G
# for lib in p.N14.plate.104720.part.1_cat p.N14.plate.104720.part.2_cat



