#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")

. ${SCRIPT_DIR}/global_config_bash.rc

if [[ $# -ne 1 ]]; then
    echo "Usage: sh ./step0b_parse_qt.sh [Manifest file Name]"
    echo "eg.: sh ./step0b_parse_qt.sh /home/wangm6/tmp2/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv"
    exit 1
else
    MANIFEST=$1
fi
DATE=`date '+%Y-%m-%d-%H-%M'`
QUEUE=seq*.q,long.q
CMD="qsub -q $QUEUE -o ${BUFFER_DIR}/logs/parse_qt_${DATE}.stdout -e ${BUFFER_DIR}/logs/parse_qt_${DATE}.stderr -N PARSE_QTs -v SCRIPT_DIR=${SCRIPT_DIR} -S /bin/sh ${SCRIPT_DIR}/parse_qt.sh $MANIFEST"
echo $CMD
eval $CMD
