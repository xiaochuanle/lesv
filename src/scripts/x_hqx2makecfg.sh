#!/bin/bash

function dumpUsage()
{
	echo "USAGE:"
	echo "$0 cfg_file"
}

if [ $# -ne 1 ];then
	dumpUsage
	exit 1
fi

CFG=$1

echo "PROJECT=" > ${CFG}
echo "RAW_READS=" >> ${CFG}
echo "REFERENCE=" >> ${CFG}
echo "TRF_FILE=" >> ${CFG}
echo "THREADS=4" >> ${CFG}

echo "" >> ${CFG}
echo "# split long read into short subreads" >> ${CFG}
echo "MAX_SUBSEQ_SIZE=50000" >> ${CFG}
echo "SUBSEQ_OVLP_SIZE=0" >> ${CFG}
echo "MIN_LAST_SUBSEQ_SIZE=20000" >> ${CFG}

echo "" >> ${CFG}
echo "# reference mapping options" >> ${CFG}
echo "MAP_OPTIONS=" >> ${CFG}

echo "" >> ${CFG}
echo "# sv read options" >> ${CFG}
echo "SVR_MIN_SEQ_SIZE=3000" >> ${CFG}
echo "SVR_MIN_SVE_PERC_IDENTITY=70.0" >> ${CFG}
echo "SVR_MAX_OVERHANG=300" >> ${CFG}

echo "" >> ${CFG}
echo "# sv signature options" >> ${CFG}
echo "SVSIG_MIN_INDEL_SIZE=40" >> ${CFG}
