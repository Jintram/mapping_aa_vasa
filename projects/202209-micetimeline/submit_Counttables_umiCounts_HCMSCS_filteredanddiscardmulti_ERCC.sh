#!/bin/sh 


################################################################################

# Note that this file overwrites the genome mapping files!

################################################################################






mappingfolder=/hpc/hub_oudenaarden/mwehrens/fastq/202209_micetimeline/mapping/

cd $mappingfolder

mem=10G
tmpspace=10G

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04/

#gtffile=/Volumes/fastq_m.wehrens/Mapping/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/filteredcustom_Homo_sapiens.GRCh38.93.gtf
    # only protein_coding, lincRNA, antisense (like cell ranger)
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107/ensembl/Mus_musculus.GRCm39.107.gtf
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCm39.107/ensembl/filteredcustom_Mus_musculus.GRCm39.107.gtf
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/ERCC/ERCC92.gtf
  # only protein_coding, lincRNA, antisense (like cell ranger)
  

#nrthreads=
trimUMI=0
#steptorun
#bamfile=..
countmulti=no

# for testing purposes:
# lib=head_p.N1.plate.97474.part.1_cat




for lib in Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_nc Martijn-Wehrens-library-test-mice-10cycles_AAC3KK3M5_S21_L001_pT Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_nc Martijn-Wehrens-library-test-mice-12cycles_AAC3KK3M5_S22_L001_pT
do
  echo $lib
  
  bamfile=${mappingfolder}/${lib}.nonRibo_E99_Aligned.out
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






