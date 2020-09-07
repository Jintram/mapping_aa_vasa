#!/bin/bash

### input paths (to modify by user)
p2s=/hpc/hub_oudenaarden/aalemany/bin/vasaplate
p2trimgalore=/hpc/hub_oudenaarden/aalemany/bin/TrimGalore-0.4.3
p2cutadapt=/hpc/hub_oudenaarden/aalemany/bin
p2bwa=/hpc/hub_oudenaarden/bin/software/bwa-0.7.10
p2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1/
p2star=/hpc/hub_oudenaarden/aalemany/bin/STAR-2.7.3a/bin/Linux_x86_64/
p2bedtools=/hpc/hub_oudenaarden/aalemany/bin/bedtools2/bin/

### check input parameters
if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) library name (prefix of the fastq files)"
    echo "2) genome: MOUSE /  HUMAN"
    echo "3) read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136)"
    echo "4) prefix for output files"
    exit
fi

lib=$1
ref=$2
n=$3
out=$4

### check existence of input fastq files
if [ ! -f ${lib}_R1.fastq.gz ]
then
    echo "${lib}_R1.fastq.gz not found"
    exit
fi

if [ ! -f ${lib}_R1.fastq.gz ]
then
    echo "${lib}_R2.fastq.gz not found"
    exit
fi

### check python version (we want version 3)
v=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))' | awk -F "." '{print $1}')
if [ $v -ne "3" ]
then
    echo "python needs to be 3"
    exit
fi

### set references
if [[ $ref == "MOUSE" ]]
then
    riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
    genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
    if [ ! -d $genome ]
    then
        echo "genome not found"
        exit
    fi
    refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
elif [[ $ref == "HUMAN" ]]
then
    riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/unique_rRNA_human.fa
    genome=/hpc/hub_oudenaarden/group_references/ensembl/99/homo_sapiens/star_v273a_NOMASK_NOERCC_index_$n
    if [ ! -d $genome ]
    then
        echo "genome not found"
        exit
    fi
    refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Homosapines_ensemble99.homemade_IntronExonTrna.bed
fi

### extract cell barcodes
jobid=extract-${lib}
jcbc=1
jcbc=$(sbatch --export=All -c 1 -N 1 -J ${jobid} -e ${jobid}.err -o ${jobid}.out -t 48:00:00 --mem=10G --wrap="${p2s}/extractBC.sh ${lib} vasaplate ${p2s}")
jcbc=$(echo $jcbc | awk '{print $NF}')
  # Some edits should be made here, eventually, the mapping should be done using 
  # both reads, and also reads from my primer could be separated to different
  # files

### trim file
jobid=trim-${lib}
jtrim=2
jtrim=$(sbatch --export=All -N 1 -J ${jobid} -e ${jobid}.err -o ${jobid}.out --dependency=afterany:$jcbc -t 15:00:00 --mem=10G --wrap="${p2s}/trim.sh ${lib}_cbc.fastq.gz ${p2trimgalore} ${p2cutadapt}")
jtrim=$(echo $jtrim | awk '{print $NF}')
  # note that this is applied to only r2 files -- i think cutadapt is compatible
  # with paired reads, but pro'lly needs extra settings

### ribo-map
jobid=ribo-${lib}
jribo=3
jribo=$(sbatch --export=All -c 1 -N 8 -J $jobid -o ${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jtrim --wrap="${p2s}/ribo-bwamem.sh $riboref ${lib}_cbc_trimmed_homoATCG.fq.gz ${lib}_cbc_trimmed_homoATCG $p2bwa $p2samtools y $p2s")
jribo=$(echo $jribo | awk '{print $NF}')
  # this is done "in silico" remove ribosomal reads (should be filtered by wet lab protocol already, but not 100%)
  # not that ribo-bwamem script also removes the reads using riboread-selection.py after mapping

### map to genome
jobid=gmap-$lib
jgmap=4
jgmap=$(sbatch --export=All -c 1 -N 8 -J $jobid -o ${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jribo --wrap="${p2s}/map_star.sh ${p2star} ${p2samtools} ${genome} ${lib}_cbc_trimmed_homoATCG.nonRibo.fastq.gz ${lib}_cbc_trimmed_homoATCG.nonRibo_E99_")
jgmap=$(echo $jgmap | awk '{print $NF}')
  # this script "just" calls the star mapping

### get bed files 
  # note that there can be ambiguities in where stuff needs to be mapped
  # mostly due to overlapping annotation for whatever reason or
  # multi-mappers
jobid=b2bs-$lib
jb2bs=5
jb2bs=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${jobid}.err -t 24:00:00 --mem=80G --dependency=afterany:$jgmap --wrap="${p2s}/deal_with_singlemappers.sh ${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y")
jb2bs=$(echo $jb2bs | awk '{print $NF}')
  # this script uses awk commands to apply a set of pre-determined rules for 
  # dealing with ambiguities in the assignment of locations to ref transcripts

jobid=b2bm-$lib
jb2bm=6
jb2bm=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jgmap --wrap="${p2s}/deal_with_multimappers.sh ${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y")
jb2bm=$(echo $jb2bm | awk '{print $NF}')
  # same, but now for multi-mappers

### count table
jobid=cout-$lib
jcout=7
jcout=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${jobid}.err -t 120:00:00 --mem=80G --dependency=afterany:$jb2bs --dependency=afterany:$jb2bm --wrap="${p2s}/countTables_final.py ${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.singlemappers_genes.bed ${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.nsorted.multimappers_genes.bed ${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned vasa")
  # then, a count table is created










