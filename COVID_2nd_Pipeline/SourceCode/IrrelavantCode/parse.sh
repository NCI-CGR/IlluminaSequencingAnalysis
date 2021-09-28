#!/bin/sh
. ${SCRIPT_DIR}/global_config_bash.rc
if [[ "$#" -ne 4 ]] && [[ "$#" -ne 3 ]]; then
    echo "Illegal number of parameters."
    echo "Usage: sh ./parse.sh /home/wangm6/tmp2/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv /home/wangm6/tmp2/AUG-FAMILIAL_manifest_update.txt redo|update [depth]"
    exit 1
else
    MANIFEST=$1
    MANIFEST_UPDATE=$2
    REDO=$3
    if [[ "$#" -eq 4 ]]; then
      DP=$4
    else
      DP=60
    fi
    if [ "$REDO" != "redo" ] && [ "$REDO" != "update" ]; then
        echo "Error: the third argument must be update or redo!"
        exit 1
    fi

 
    DATE=`date +%Y%m%d`
    OUTPUT_DIR=${BUFFER_DIR}/output/$DATE
fi


echo $OUTPUT_DIR
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p $OUTPUT_DIR
fi
# CMD="java -jar secondaryParsing.jar parse $MANIFEST ${SCRIPT_DIR}/input/all_merging_records.txt ${SCRIPT_DIR}/step2_merge.sh ${OUTPUT_DIR}/manifest_errs.txt ${OUTPUT_DIR}/new_update_analysisIDs.txt ${OUTPUT_DIR}/filenames_restore.txt 60 update dedup"

CMD="java -jar ${SCRIPT_DIR}/secondaryParsing.jar parse $MANIFEST ${BUFFER_DIR}/input/all_merging_records.txt ${OUTPUT_DIR}/step2_merge_${DATE}.sh ${OUTPUT_DIR}/manifest_errs.txt ${OUTPUT_DIR}/new_update_analysisIDs.txt ${OUTPUT_DIR}/filenames_restore.txt ${MANIFEST_UPDATE} $DP $REDO ${BUFFER_DIR}/input/qt_all.txt"

# CMD="java -jar secondaryParsing.jar parse $MANIFEST ${SCRIPT_DIR}/input/all_merging_records.txt ${SCRIPT_DIR}/step2_merge_${DATE}.sh ${OUTPUT_DIR}/manifest_errs.txt ${OUTPUT_DIR}/new_update_analysisIDs.txt ${OUTPUT_DIR}/filenames_restore.txt ${MANIFEST_UPDATE} 60 update paired"
echo $CMD
eval $CMD
