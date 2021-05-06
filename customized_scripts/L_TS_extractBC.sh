#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give: "
    echo "1) path to general parameter file"
    echo "2) path to run parameter file"
    echo "3) protocol [celseq1, celseq2, vasaplate]"
    echo "4) library prefix"
    exit
fi

general_parameter_filepath=$1
run_parameter_filepath=$2
protocol=$3
lib=$4

source $general_parameter_filepath
source $run_parameter_filepath

# Todo: below, put "-pf ${path2scripts}/targetedprimers.tsv" in parameter that's dependent on whether we do TS or not
if [ $TS = '1' ]; then
    TS_str="-TS 1 -pf ${path2scripts}/targetedprimers.tsv"
else
    TS_str=""
fi

if [ $protocol == 'celseq1' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq1.tsv --cbchd 0 --lenumi 4 --outdir ${outdir} $TS_str
elif [ $protocol == 'celseq2' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --outdir ${outdir} $TS_str
elif [ $protocol == 'vasaplate' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --outdir ${outdir} $TS_str
else
    echo "unknown protocol [celseq1, celseq2, vasaplate]"
    exit
fi
