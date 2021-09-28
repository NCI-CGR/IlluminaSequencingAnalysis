#!/bin/sh
. ${SCRIPT_DIR}/global_config_bash.rc

if [[ $# -ne 1 ]]; then
    echo "Usage: sh ./parse_qt.sh [Manifest file Name]"
    echo "eg.: sh ./parse_qt.sh /home/wangm6/tmp2/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv"
    exit 1
else
    MANIFEST=$1
fi
DATE=`date '+%Y%m%d'`
if [[ ! -d  ${BUFFER_DIR}/output/${DATE} ]]; then
   mkdir -p ${BUFFER_DIR}/output/${DATE}
   if [[ $? -ne 0 ]]; then
      echo "Error: failed to create the folder ${BUFFER_DIR}/output/${DATE}!"
      exit 1
   fi
fi
CMD="java -jar ${SCRIPT_DIR}/secondaryParsing.jar qt $MANIFEST ${BUFFER_DIR}/input/qt_all.txt  ${BUFFER_DIR}/output/${DATE}/fastq_restoration.txt  ${BUFFER_DIR}/output/${DATE}/samples_retrimming.txt ${BUFFER_DIR}/output/${DATE}/qt_err.txt ${BUFFER_DIR}/output/${DATE}/samples_new_trimmed.txt"
echo $CMD
eval $CMD
