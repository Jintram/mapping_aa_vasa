#!/bin/sh 

cd /hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping.93.may25/

mem=10G
tmpspace=10G

script_path=/hpc/hub_oudenaarden/mwehrens/scripts/mapping_2021-04/

#gtffile=/Volumes/fastq_m.wehrens/Mapping/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
#gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/Homo_sapiens.GRCh38.93.gtf
gtffile=/hpc/hub_oudenaarden/mwehrens/ref/GRCh38.93/ensembl/filteredcustom_Homo_sapiens.GRCh38.93.gtf
  # only protein_coding, lincRNA, antisense (like cell ranger)

#nrthreads=
trimUMI=0
#steptorun
#bamfile=..
countmulti=no

# for testing purposes:
# lib=head_p.N1.plate.97474.part.1_cat

for lib in JE7_AHFL7NBGX5_S16_cat_nc JE7_AHFL7NBGX5_S16_cat_pT JE6_AHFL77BGX5_S7_cat_nc JE6_AHFL77BGX5_S7_cat_pT HUB-AL-s002_HG25TBGXF_S6_cat_nc HUB-AL-s002_HG25TBGXF_S6_cat_pT JE8_AHFL7NBGX5_S17_cat_nc JE8_AHFL7NBGX5_S17_cat_pT HUB-AL-s001_HG25TBGXF_S5_cat_nc HUB-AL-s001_HG25TBGXF_S5_cat_pT JE4_AHFL7NBGX5_S4_cat_nc JE4_AHFL7NBGX5_S4_cat_pT HUB-MW-006_AH32W2BGX9_S6_cat_nc HUB-MW-006_AH32W2BGX9_S6_cat_pT HUB-MW-008_HC3GFBGX9_S7_cat_nc HUB-MW-008_HC3GFBGX9_S7_cat_pT HUB-MW-005_AH32W2BGX9_S5_cat_nc HUB-MW-005_AH32W2BGX9_S5_cat_pT HUB-MW-007_HC3GFBGX9_S6_cat_nc HUB-MW-007_HC3GFBGX9_S6_cat_pT rJE1_AHFL7NBGX5_S3_cat_nc rJE1_AHFL7NBGX5_S3_cat_pT HUB-JE-011_HGVN3BGX9_S2_cat_nc HUB-JE-011_HGVN3BGX9_S2_cat_pT JE5_AHFL77BGX5_S6_cat_nc JE5_AHFL77BGX5_S6_cat_pT JE3_AHY3WGBGX3_S2_cat_nc JE3_AHY3WGBGX3_S2_cat_pT JE2_AHY3WGBGX3_S1_cat_nc JE2_AHY3WGBGX3_S1_cat_pT HUB-JE-010_HGVN3BGX9_S1_cat_nc HUB-JE-010_HGVN3BGX9_S1_cat_pT head_HUB-JE-010_HGVN3BGX9_S1_cat_nc head_HUB-JE-010_HGVN3BGX9_S1_cat_pT 
do
  echo $lib
  
  bamfile=/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping.93.may25/${lib}.nonRibo_E99_Aligned.out
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






