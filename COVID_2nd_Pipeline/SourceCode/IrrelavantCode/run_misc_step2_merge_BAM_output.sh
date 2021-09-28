#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")

. ${SCRIPT_DIR}/global_config_bash.rc
if [[ $# -ne 1 ]]; then
  echo "Error: please input build folder name"
  exit 1
fi
BUILD_DIR=$1
cd ${BUILD_DIR}
BUILD_DIR=`pwd`
echo $BUILD_DIR
BUILD_NAME=`echo ${BUILD_DIR} | rev | cut -d"/" -f 1 | rev`
echo $BUILD_NAME
LOG_DIR=${BUFFER_DIR}/logs/parse_${BUILD_NAME}
OUT_DIR=${BUFFER_DIR}/output/bam_results_${BUILD_NAME}


CMD="${SCRIPT_DIR}/misc_merge_all.sh ${OUT_DIR}"
SGE_CMD="qsub -q all.q -cwd -j y -o ${BUFFER_DIR}/logs/merge_bam_output_${BUILD_NAME}.stdout -N merge_output -v SCRIPT_DIR=$SCRIPT_DIR -S /bin/sh $CMD"
echo $SGE_CMD
eval $SGE_CMD
echo "Please check the results at: ${OUT_DIR}/report_cmp.txt ${OUT_DIR}/err_cmp.txt"
