#!/bin/sh

#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")
. ${SCRIPT_DIR}/global_config_bash.rc
QUEUE=seq*.q,long.q
DATE=`date '+%Y-%m-%d-%H-%M'`

rm -f ${BUFFER_DIR}/logs/scan_merging_records_${DATE}.*
CMD="qsub -q $QUEUE -o ${BUFFER_DIR}/logs/scan_merging_records_${DATE}.stdout -e ${BUFFER_DIR}/logs/scan_merging_records_${DATE}.stderr -N SCAN_LOGs -v SCRIPT_DIR=$SCRIPT_DIR -S /bin/sh ${SCRIPT_DIR}/scan_logs.sh"
echo $CMD
eval $CMD
