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

CONFIG=$1

while read line;do  
    eval "$line"  
done < ${CONFIG}

CNS_SVR_GROUP_DIR="${PROJECT}/6-svr-group-consensus"

if [ ! -d ${CNS_SVR_GROUP_DIR} ]; then
	echo "Directory '${CNS_SVR_GROUP_DIR}' does not exist"
	exit 1
fi

SV_CALL_DIR="${PROJECT}/7-sv-call"
if [ ! -d ${SV_CALL_DIR} ]; then
	mkdir -p ${SV_CALL_DIR}
fi
SV_CALL_DONE="${SV_CALL_DIR}/sv_call.done"
if [ -f ${SV_CALL_DONE} ]; then
	echo "[${SV_CALL_DONE}] is done."
	exit 0
fi

BAM_DIR="${SV_CALL_DIR}/bam"
if [ ! -d ${BAM_DIR} ]; then
	mkdir -p ${BAM_DIR}
fi

SORTED_BAM_DIR="${SV_CALL_DIR}/sorted_bam"
if [ ! -d ${SORTED_BAM_DIR} ]; then
	mkdir -p ${SORTED_BAM_DIR}
fi

SORT_BAM_DONE="${BAM_DIR}/view_and_sort_bam.done"
if [ -f ${SORT_BAM_DONE} ]; then
	echo "[view and sort bam] is done."
else
	sam_file_extension=(".sam")
	for file in ` ls ${CNS_SVR_GROUP_DIR} `
	do
		fileName=${CNS_SVR_GROUP_DIR}"/"$file
		EXTENSION="."${fileName##*.}
		if [[ "${sam_file_extension[@]/$EXTENSION/}" != "${sam_file_extension[@]}" ]];then
			echo "Process ${fileName}"
			CMD="samtools view -bS ${fileName}"
			bam_file="${BAM_DIR}/${file}.bam"
			${CMD} > ${bam_file}
			if [ $? -ne 0 ]; then
				echo "[${CMD}] Fail"
				exit 1
			fi
			sorted_bam_file_prefix="${SORTED_BAM_DIR}/${file}"
			CMD="samtools sort -@${THREADS} ${bam_file} ${sorted_bam_file_prefix}"
			${CMD}
			if [ $? -ne 0 ]; then
				echo "[${CMD}] Fail."
				exit 1
			fi
		fi
	done
	touch ${SORT_BAM_DONE}
fi

PBSV_DIR="${SV_CALL_DIR}/pbsv"
if [ ! -d ${PBSV_DIR} ]; then
	mkdir -p ${PBSV_DIR}
fi
PBSV_DISCOVER_DONE="${PBSV_DIR}/pbsv_discover.done"
if [  -f ${PBSV_DISCOVER_DONE} ]; then
	echo "[pbsv signature discover] is done."
else
        bam_file_extension=(".bam")
        for file in ` ls ${SORTED_BAM_DIR} `
        do
                fileName=${SORTED_BAM_DIR}"/"$file
                EXTENSION="."${fileName##*.}
                if [[ "${bam_file_extension[@]/$EXTENSION/}" != "${bam_file_extension[@]}" ]];then
			CMD="pbsv discover ${fileName} ${fileName}.svsig.gz -l 20"
			if [ -n "${TRF_FILE}" ]; then
				CMD="${CMD} --tandem-repeats ${TRF_FILE}"
			fi
			echo "signature discover for ${file}"
			${CMD}
			if [ $? -ne 0 ]; then
				echo "[${CMD}] Fail."
				exit 1
			fi
                fi
        done	
	touch ${PBSV_DISCOVER_DONE}
fi

PBSV_CALL_DONE="${PBSV_DIR}/pbsv_call.done"
if [ -f ${PBSV_CALL_DONE} ]; then
	echo "[pbsv call sv] is done."
else
	CMD="pbsv call ${REFERENCE} ${SORTED_BAM_DIR}/*.svsig.gz ${PBSV_DIR}/${PROJECT}.pbsv.vcf -t INS,DEL --max-ins-length 30k"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	CMD="bgzip ${PBSV_DIR}/${PROJECT}.pbsv.vcf"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	CMD="tabix ${PBSV_DIR}/${PROJECT}.pbsv.vcf.gz"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[${CMD}] Fail."
		exit 1
	fi
	touch ${PBSV_CALL_DONE}
fi
