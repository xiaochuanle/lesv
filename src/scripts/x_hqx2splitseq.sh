#!/bin/bash

function dumpUsage()
{
	echo "USAGE:"
	echo "$0 work_directory max_subseq_length subseq_ovlp_length min_last_subseq_length input output"
}

if [ $# -ne 6 ];then
	dumpUsage
	exit 1
fi

WRK_DIR=$1
MAX_L=$2
OVLP_L=$3
LAST_L=$4
INPUT=$5
OUTPUT=$6

if [ ! -f ${WRK_DIR} ]; then
	mkdir -p ${WRK_DIR}
fi

JOB_DONE="${WRK_DIR}/split_reads.done"

if [ -f ${JOB_DONE} ]; then
	echo "[$0] has been done."
	exit 0
fi

CMD="qx2splitseq ${MAX_L} ${OVLP_L} ${LAST_L} ${INPUT}"

$CMD > ${OUTPUT}

if [ $? -ne 0 ]; then
	echo "[$CMD] Fail"
	exit 1
fi

touch ${JOB_DONE}
exit 0
