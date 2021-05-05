#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give: "
    echo "1) path to run parameter file"
    echo "2) protocol [celseq1, celseq2, vasaplate]"
    exit
fi

run_parameter_filepath=$1
protocol=$2

source ./set_paths.sh
source $run_parameter_filepath

if [ $protocol == 'celseq1' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq1.tsv --cbchd 0 --lenumi 4 --outdir ${outdir} -pf ${path2scripts}/targetedprimers.tsv
elif [ $protocol == 'celseq2' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --outdir ${outdir} -pf ${path2scripts}/targetedprimers.tsv
elif [ $protocol == 'vasaplate' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --outdir ${outdir} -pf ${path2scripts}/targetedprimers.tsv
else
    echo "unknown protocol [celseq1, celseq2, vasaplate]"
    exit
fi
