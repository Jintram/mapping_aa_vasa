#!/bin/bash

email=a.alemany@hubrecht.eu
Dt=15:00:00
Dmem=10G
threads=1

if [ $# -ne 2 ]
then
    echo "Please, give these following inputs:"
    echo "1) input root to fastq files"
    echo "2) protocol [celseq1, celseq2, scscar, nla]"
    exit
fi

infq=$1
protocol=$2

source /hpc/hub_oudenaarden/aalemany/virtualEnvironments/venv36/bin/activate
echo "/hpc/hub_oudenaarden/aalemany/bin/mapandgo2/extractBC.sh $infq $protocol" | qsub -V -cwd -N bc-${infq} -o bc-${infq}.out -e bc-${infq}.err -m eas -M ${email} -pe threaded ${threads} -l h_rt=${Dt} -l h_vmem=${Dmem}


