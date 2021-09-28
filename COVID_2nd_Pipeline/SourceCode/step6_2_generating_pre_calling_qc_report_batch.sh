#!/bin/sh

. ./global_config_bash.rc

#find the most recent coverage report file as the coverage report file. There is a potential bug that it only works on the most recently created report file. 
#COVERAGE_REPORT_FILE=`find ${COVERAGE_REPORT_DIR} -name coverage_report_????????.txt  -print0 | xargs -0r ls -ltr | tail -n 1 | awk '{print $9}'`
PREQC_REPORT_FILE=$1
DATE=$(basename $PREQC_REPORT_FILE .txt | cut -f5 -d_)

LOG_DIR=${CLUSTER_JOB_LOG_DIR}/${DATE}

for SINGLE_OUT in ${LOG_DIR}/_pre_calling_qc_report_*.stdout;do
#	REPORT_LINE=$(grep -A 1 "Total Exome Bases >= 50x" $SINGLE_OUT | tail -n 1)
#	echo -e $REPORT_LINE  >> $COVERAGE_REPORT_FILE
	grep -A 1 "Contents to be streamed to the flowcell-level report file" $SINGLE_OUT | tail -n 1 >> $PREQC_REPORT_FILE
done

