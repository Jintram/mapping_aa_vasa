#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give: "
    echo "1) path to general parameter file"
    echo "2) path to run parameter file"
    echo "3) protocol [celseq1, celseq2, vasaplate, takara]"
    echo "4) library prefix"
    exit
fi

general_parameter_filepath=$1
run_parameter_filepath=$2
protocol=$3
lib=$4

source $general_parameter_filepath
source $run_parameter_filepath

if [[ $TMPDIR == "" ]]; then
  TMPDIR=$outdir
fi

# Todo: below, put "-pf ${path2scripts}/targetedprimers.tsv" in parameter that's dependent on whether we do TS or not
if [ $TS = '1' ]; then
    TS_str="-TS 1 -pf ${path2scripts}/targetedprimers.tsv"
else
    TS_str=""
fi

if [ $protocol == 'celseq1' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq1.tsv --cbchd 0 --lenumi 4 --outdir ${TMPDIR} $TS_str
    exitcode=$?
elif [ $protocol == 'celseq2' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --TTTfilter --outdir ${TMPDIR} $TS_str
    exitcode=$?
elif [ $protocol == 'vasaplate' ]
then
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --TTTfilter --outdir ${TMPDIR} $TS_str
    $exitcode = $?
elif [ $protocol == 'takara' ]
then
    # note that this a bit of boiler plating, tailored for the Wang et al. paper
    # (Wang, L., Yu, P., Zhou, B., Song, J., Li, Z., Zhang, M., … Hu, S. (2020). Single-cell reconstruction of the adult human heart during heart failure and recovery reveals the cellular landscape underlying cardiac function. Nature Cell Biology. https://doi.org/10.1038/s41556-019-0446-7)
    # They use icell8 platform, and when I look at their processed and meta data (GSE109816, GSE121893), it appears
    # a barcode of 11, and a UMI of 11-14 [I think this is 14, but limited by R1 seq length]
    ${pythonbin} ${path2scripts}/TS_concatenator.py --fqf ${lib} --cbcfile ${path2scripts}/bc_takara_map_tcr_v1.tsv --cbchd 0 --lenumi 14 --lencbc 11 --outdir ${TMPDIR} $TS_str
    exitcode=$?
      # Note: by *not* giving the option "--umifirst", it will be set to 0.
      # idem for  --TTTfilter (which is non-consequential, given that this data won't contain polyT seqs anyhow.)
else
    echo "unknown protocol [celseq1, celseq2, vasaplate]"
    exit
fi

# let's copy some files to output dir in case we're working in temporary directory
if [[ ${TMPDIR} != ${outdir} ]]; then
  if [[ -f ${TMPDIR}/${lib}_pT_R2_cbc.fastq.gz ]]; then
    cp ${TMPDIR}/${lib}_pT_R2_cbc.fastq.gz $outdir
  else
    echo "${TMPDIR}/${lib}_pT_R2_cbc.fastq.gz not found"
  fi
  if [[ -f ${TMPDIR}/${lib}_nc_R2_cbc.fastq.gz ]]; then
    cp ${TMPDIR}/${lib}_nc_R2_cbc.fastq.gz $outdir
  else
    echo "${TMPDIR}/${lib}_nc_R2_cbc.fastq.gz not found"
  fi
  echo "(todo) _TS file will be added later"
fi

# if the script didn't end normally, also throw error here
if [[ $exitcode -ne 0 ]]; then
  echo "Non-proper ending of script detected."
  exit 1
fi
















