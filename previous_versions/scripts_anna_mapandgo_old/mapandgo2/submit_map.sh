#!/bin/bash

email=a.alemany@hubrecht.eu
Dt=05:00:00
Dmem=40G
threads=6

if [ $# -ne 3 ]
then
    echo "Please, give"
    echo "1) root input fastq lanes"
    echo "2) protocol"
    echo "3) reference"
    exit
fi

echo "/hpc/hub_oudenaarden/aalemany/bin/mapandgo2/map.sh $1 $2 $3" | qsub -cwd -N map-$1 -o map-${1}.out -e map-${1}.err -m eas -M ${email} -pe threaded ${threads} -l h_rt=${Dt} -l h_vmem=${Dmem} -l hostname='!n008*'




