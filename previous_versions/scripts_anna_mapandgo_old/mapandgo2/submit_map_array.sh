#!/bin/bash

p2s=/hpc/hub_oudenaarden/aalemany/bin/mapandgo2
email=a.alemany@hubrecht.eu

if [ $# -ne 1 ]
then
    echo "Please, give input root file"
    exit
fi

in=$1
out=${in%_*_S*_L*}

# merge
Dt=15:00:00
Dmem=10G
name=merge-${out}
threads=2
echo "${p2s}/mergeLanes.sh $1 $2" | qsub -cwd -N $name -o ${name}.out -e ${name}.err -m eas -M ${email} -pe threaded ${threads} -l h_rt=${Dt} -l h_vmem=${Dmem}

# extract barcodes
name=extract-$out
echo "${p2s}/extractBC.sh $out celseq2" | qsub -hold_jid merge-${out} -cwd -N $name -o ${name}.out -e ${name}.err -m eas -M ${email} -pe threaded ${threads} -l h_rt=${Dt} -l h_vmem=${Dmem}
