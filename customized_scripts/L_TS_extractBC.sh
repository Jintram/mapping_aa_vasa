#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give: "
    echo "1) input root to fastq files"
    echo "2) protocol [celseq1, celseq2, vasaplate]"
    echo "3) path to concatenator.py"
    exit
fi

outfq=$1
protocol=$2
path2scripts=$3
outdir=./
pythonbin=/Users/m.wehrens/anaconda3/bin/python

if [ $protocol == 'celseq1' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq1.tsv --cbchd 0 --lenumi 4 --outdir ${outdir} -pf ${path2scripts}/targetedprimers.tsv
elif [ $protocol == 'celseq2' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --outdir ${outdir} -pf ${path2scripts}/targetedprimers.tsv
elif [ $protocol == 'vasaplate' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --outdir ${outdir} -pf ${path2scripts}/targetedprimers.tsv
else
    echo "unknown protocol [celseq1, celseq2, vasaplate]"
    exit
fi
