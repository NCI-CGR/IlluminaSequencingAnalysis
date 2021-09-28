#!/bin/sh

# To scan all Quality trimming logs after a specific date (or all logs) to save the results
# E.g. command line: 
# sh ./step0_run_scan_quality_trimming.sh 2017-12-01

SCRIPT=$(readlink -f "$0")
SCRIPT_HOME=$(dirname "$SCRIPT")

. ${SCRIPT_HOME}/global_config_bash.rc
if [[ ! -d ${BUFFER_DIR}/logs ]]; then
  mkdir -p ${BUFFER_DIR}/logs
fi 
if [[ $# -eq 1 ]]; then
  DATE_NEWER=$1
else
  if [[ $# -gt 1 ]]; then
    echo "Usage: sh ./step0_run_scan_quality_trimming.sh [date]"
    echo "Scan all quality trimming(QT) logs which are newer than speicfied date in primary analysis to get all the new QT list"
    echo "E.g., sh ./step0_run_scan_quality_trimming.sh 2016-17-01"
    exit 1
  fi
fi
QUEUE=seq*.q
DATE=`date '+%Y-%m-%d-%H-%M'`
rm -f ${SCRIPT_HOME}/logs/scan_qt_${DATE}.std???
if [[ $# -eq 1 ]]; then
  CMD="qsub -q $QUEUE -o ${BUFFER_DIR}/logs/scan_qt_${DATE}.stdout -e ${BUFFER_DIR}/logs/scan_qt_${DATE}.stderr -N SCAN_QTs -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/scan_quality_trimming.sh $DATE_NEWER"
else
  CMD="qsub -q $QUEUE -o ${BUFFER_DIR}/logs/scan_qt_${DATE}.stdout -e ${BUFFER_DIR}/logs/scan_qt_${DATE}.stderr -N SCAN_QTs -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/scan_quality_trimming.sh"
fi
echo $CMD
eval $CMD
