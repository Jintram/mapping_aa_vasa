#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give input (1) root file; (2) protocol [celseq1, celseq2, nla, scarsc]; (3) reference [mouse, human, zebrafish, spiny]"
    exit
fi

p2s=/hpc/hub_oudenaarden/aalemany/bin/mapandgo2
p2trimgalore=/hpc/hub_oudenaarden/aalemany/bin/TrimGalore-0.4.3
p2cutadapt=/hpc/hub_oudenaarden/aalemany/bin
p2star=/hpc/hub_oudenaarden/avo/nascent/STAR-2.5.3a/bin/Linux_x86_64
p2bedtools=/hpc/hub_oudenaarden/aalemany/bin/bedtools2/bin
p2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1

source /hpc/hub_oudenaarden/aalemany/virtualEnvironments/venv36/bin/activate

in=$1
out=${in%_*_S*_L*}
protocol=$2
reference=$3

email=a.alemany@hubrecht.eu

# merge data
echo "${p2s}/mergeLanes.sh $in $out" | qsub -cwd -m eas -M $email -N merge-$out -e merge-${out}.err -o merge-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -pe threaded 2

# extract barcodes
echo "${p2s}/extractBC.sh $out ${protocol} ${p2s}" | qsub -V -cwd -m eas -M $email -N extract-$out -e extract-${out}.err -o extract-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -hold_jid merge-$out

# trim
echo "${p2s}/trim.sh ${out}_cbc.fastq.gz ${p2trimgalore} ${p2cutadapt}" | qsub -cwd -m eas -M $email -N trim-$out -e trim-${out}.err -o trim-${out}.out -l h_vmem=10G -l h_rt=15:00:00 -hold_jid extract-$out

# map with star
echo "${p2s}/mapstar.sh ${out}_cbc_trimmed.fq.gz ${out}_cbc_trimmed_star $reference ${p2star}" | qsub -cwd -m eas -M $email -N map-$out -e map-${out}.err -o map-${out}.out -l h_vmem=30G -l h_rt=15:00:00 -pe threaded 12 -hold_jid trim-$out

# create count tables from star map
echo "${p2s}/getIntronsExons.sh ${out}_cbc_trimmed_starAligned.sortedByCoord.out.bam $reference ${out}_cbc_trimmed_star ${p2bedtools} ${p2samtools} ${p2s}" | qsub -V -cwd -m eas -M $email -N tab-$out -e tab-${out}.err -o tab-${out}.out -l h_vmem=20G -l h_rt=15:00:00 -pe threaded 2 -hold_jid map-$out




