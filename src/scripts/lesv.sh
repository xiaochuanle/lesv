#!/bin/bash

function dumpUsage()
{
	echo "USAGE:"
	echo "$0 config|run cfg_file"
}

if [ $# -ne 2 ];then
	dumpUsage
	exit 1
fi

TASK=$1
CONFIG=$2

if [ "$TASK" == "config" ]; then
	x_hqx2makecfg.sh $CONFIG
	exit 0
elif [ "$TASK" != "run" ]; then
	dumpUsage
	exit 1
fi
	
	

while read line;do  
    eval "$line"  
done < ${CONFIG}  

function validateArguments()
{
	if [ ! -n "${PROJECT}" ]; then
		echo "Argument to option [PROJECT] in config file [${CONFIG}] must be specified"
		exit 1
	fi

	if [ ! -n "${RAW_READS}" ]; then
		echo "Argument to option [RAW_READS] in config file [${CONFIG}] must be specified"
		exit 1
	fi

	if [ ! -n "${REFERENCE}" ]; then
		echo "Argument to option [REFERENCE] in config file [${CONFIG}] must be specified"
		exit 1
	fi

	if [ ! -n "${THREADS}" ]; then
		echo "[WARNING] Argument to option [THREADS] in config file [${CONFIG}] is not specified, we set it to 1"
		THREADS=1
	fi

	if [ ! -n "${MAX_SUBSEQ_SIZE}" ]; then
		echo "Argument to option [MAX_SUBSEQ_SIZE] in config file [${CONFIG}] must be specified"
		exit 1
	fi

	if [ ! -n "${SUBSEQ_OVLP_SIZE}" ]; then
		echo "Argument to option [SUBSEQ_OVLP_SIZE] in config file [${CONFIG}] must be specified"
		exit 1
	fi

	if [ ! -n "${MIN_LAST_SUBSEQ_SIZE}" ]; then
		echo "Argument to option [MIN_LAST_SUBSEQ_SIZE] in config file [${CONFIG}] must be specified"
		exit 1
	fi
}

validateArguments

LAST_DIR=$(pwd)
echo ${LAST_DIR}

if [ ! -d ${PROJECT} ]; then
	mkdir -p ${PROJECT}
fi

### step 1: split reads

RAW_DATA_DIR=${PROJECT}/1-raw-data
if [ ! -d ${RAW_DATA_DIR} ]; then
	mkdir -p ${RAW_DATA_DIR}
fi

SPLIT_READS="${RAW_DATA_DIR}/split_reads.fasta"
job_worker=x_hqx2splitseq.sh
CMD="${job_worker} ${RAW_DATA_DIR} ${MAX_SUBSEQ_SIZE} ${SUBSEQ_OVLP_SIZE} ${MIN_LAST_SUBSEQ_SIZE} ${RAW_READS} ${SPLIT_READS}"
$CMD
if [ $? -ne 0 ]; then
	echo "[${CMD}] Fail"
	exit 1
fi

### step 2: reference mapping

MAP_DIR=${PROJECT}/2-map
if [ ! -d ${MAP_DIR} ]; then
	mkdir -p ${MAP_DIR}
fi

PDB_DIR=${MAP_DIR}/packed_db
RM_RESULTS=${MAP_DIR}/map.m4
job_worker=qx2map
CMD="${job_worker} ${MAP_OPTIONS} -db_dir ${PDB_DIR} -keep_db -num_threads ${THREADS} ${SPLIT_READS} ${REFERENCE}"
RM_DONE="${MAP_DIR}/map.done"
if [ -f ${RM_DONE} ]; then
	echo "[${job_worker}] is done."
else
	$CMD > ${RM_RESULTS}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail"
		exit 1
	fi
	touch ${RM_DONE}
fi

### transfer m4 to binary

M4X_DONE="${PDB_DIR}/m4x.done"
job_worker=qx2m4x
if [ -f ${M4X_DONE} ]; then
	echo "[${job_worker}] is done."
else
	CMD="${job_worker} ${PDB_DIR}"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail"
		exit 1
	fi
	touch ${M4X_DONE}
fi

### step 3: find sv reads

SVR_DIR="${PROJECT}/3-sv-reads"
if [ ! -d ${SVR_DIR} ]; then
	mkdir -p ${SVR_DIR}
fi

job_worker=qx2svr
SVR_DONE="${SVR_DIR}/svr.done"
if [ -f ${SVR_DONE} ]; then
	echo "[${job_worker}] is done."
else
	CMD="${job_worker} ${SVR_DIR} ${PDB_DIR} ${SVR_MIN_SEQ_SIZE} ${SVR_MIN_SVE_PERC_IDENTITY} ${SVR_MAX_OVERHANG} ${THREADS} ${TRF_FILE}"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	touch ${SVR_DONE}
fi

### step 4: find sv signature

SVSIG_DIR="${PROJECT}/4-sv-signature"
if [ ! -d ${SVSIG_DIR} ]; then
	mkdir -p ${SVSIG_DIR}
fi
SVSIG_DONE="${SVSIG_DIR}/svsig.done"
job_worker=qx2svsig
if [ -f ${SVSIG_DONE} ]; then
	echo "[${job_worker}] is done."
else
	CMD="${job_worker} ${SVSIG_DIR} ${SVR_DIR} ${PDB_DIR} ${SVSIG_MIN_INDEL_SIZE} ${THREADS} ${TRF_FILE}"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	touch ${SVSIG_DONE}
fi

### step 5: make sv read groups
SVR_GROUP_DIR="${PROJECT}/5-svr-group"
if [ ! -d ${SVR_GROUP_DIR} ]; then
	mkdir -p ${SVR_GROUP_DIR}
fi
SVR_GROUP_DONE="${SVR_GROUP_DIR}/svr_group.done"
job_worker=qx2msvrg
if [ -f ${SVR_GROUP_DONE} ]; then
	echo "[${job_worker}] is done."
else
	CMD="${job_worker} ${SVR_GROUP_DIR} ${SVSIG_DIR} ${PDB_DIR}"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	touch ${SVR_GROUP_DONE}
fi

### step 6: correct sv read grops
SVR_GROUP_CNS_DIR="${PROJECT}/6-svr-group-consensus"
if [ ! -d ${SVR_GROUP_CNS_DIR} ]; then
	mkdir -p ${SVR_GROUP_CNS_DIR}
fi
SVR_GROUP_CNS_DONE="${SVR_GROUP_CNS_DIR}/svr_group_cns.done"
job_worker=qx2csvrg
if [ -f ${SVR_GROUP_CNS_DONE} ]; then
	echo "[${job_worker}] is done."
else
	CMD="${job_worker} -num_threads ${THREADS} ${SVR_GROUP_CNS_DIR} ${SVR_GROUP_DIR} ${PDB_DIR}"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	touch ${SVR_GROUP_CNS_DONE}
fi

### step 7: map corrected sv reads
CNS_SVR_MAP_DONE="${SVR_GROUP_CNS_DIR}/cns_svr_mapped.done"
job_worker=qx2asvr
if [ -f ${CNS_SVR_MAP_DONE} ]; then
	echo "[${job_worker}] is done."
else
	CMD="${job_worker} ${PDB_DIR} ${SVR_GROUP_CNS_DIR} ${THREADS}"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	touch ${CNS_SVR_MAP_DONE}
fi

## step 7: call sv

CMD="x_hqx2callsv.sh ${CONFIG}"
echo ${CMD}
${CMD}
if [ $? -ne 0 ]; then
	echo "[${CMD}] Fail."
	exit 1
fi
