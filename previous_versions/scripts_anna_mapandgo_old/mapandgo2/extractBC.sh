#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give: "
    echo "1) input root to fastq files"
    echo "2) protocol [celseq1, celseq2, scscar, nla]"
    echo "3) path to concatenator.py"
    exit
fi

outfq=$1
protocol=$2
path2scripts=$3

if [ $protocol == 'celseq1' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq1.tsv --cbchd 0 --lenumi 4
    gzip ${outfq}_cbc.fastq
elif [ $protocol == 'celseq2' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst
    gzip ${outfq}_cbc.fastq
elif [ $protocol == 'scscar' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_scarsc.tsv --cbchd 0 --lenumi 3 --umifirst
    gzip ${outfq}_cbc.fastq
elif [ $protocol == 'nla' ]
then
    python3 ${path2scripts}/concatenator.py --fqf ${outfq} --cbcfile ${path2scripts}/bc_nla.tsv --cbchd 0 --lenumi 3 --umifirst
    gzip ${outfq}_cbc.fastq
else
    echo "unknown protocol [celseq1, celseq2, scscar, nla]"
    exit
fi
