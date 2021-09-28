#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")

. ${SCRIPT_DIR}/global_config_bash.rc
QUEUE=long.q
if [[ $# -ne 2 ]] && [[ $# -ne 3 ]]; then
   echo "Usage: sh ./step1_run_parse.sh [MANIFEST csv file] [update|redo] [depth]"
   echo "Eg. sh ./step1_run_parse.sh /home/wangm6/tmp2/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv update 60"
   exit 1
else
   MANIFEST=$1
   SCRIPT_MODE=$2
   if [[ $# -eq 3 ]]; then
     DP=$3
   else
     DP=60
   fi
fi
DIRNAME=`dirname $MANIFEST`
BASENAME=`basename $MANIFEST .csv`
DATE=`date '+%Y-%m-%d-%H-%M'`
rm -f ${BUFFER_DIR}/logs/step1_parse${DATE}.*
CMD="qsub -q $QUEUE -o ${BUFFER_DIR}/logs/step1_parse${DATE}.stdout -e ${BUFFER_DIR}/logs/step1_parse${DATE}.stderr -N parse_manifest -v SCRIPT_DIR=$SCRIPT_DIR -S /bin/sh ${SCRIPT_DIR}/parse.sh $MANIFEST ${DIRNAME}/${BASENAME}_update.txt $SCRIPT_MODE $DP"
echo $CMD
eval $CMD
