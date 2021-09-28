#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")
. ${SCRIPT_DIR}/global_config_bash.rc

if [[ $# -ne 1 ]] ; then 
    echo "Error: please input bam_location folder name and build name!"
    echo "E.g., /DCEG/Projects/Exome/builds/build_SR0403-001_Mexican_breast_cancer_2019_22639"
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
COUNT=`ls ${BUILD_DIR}/Manifest/*.csv | wc -l`
if [[ $COUNT -eq 1 ]]; then
  mkdir -p $OUT_DIR
  mkdir -p $LOG_DIR
  MANIFEST=`ls ${BUILD_DIR}/Manifest/*.csv`
  rm -f ${OUT_DIR}/manifest.csv
  cp -l $MANIFEST ${OUT_DIR}/manifest.csv
  if [[ $? -ne 0 ]]; then
    echo "Error: failed in cp -l $MANIFEST ${OUT_DIR}/manifest.csv"
    exit 1
  fi
else
  echo "Error: there are more than one manifest or no manifest in ${BUILD_DIR}/Manifest"
  exit 1
fi
cd ${BUILD_DIR}/bam_location
for i in */*.bam; do
  BASE_NAME=`basename $i .bam`
  OUTPUT_FILE=${OUT_DIR}/${BASE_NAME}.out
  CMD="${SCRIPT_DIR}/misc_bam_check2.sh ${BUILD_DIR}/bam_location/$i $OUTPUT_FILE $BASE_NAME"
  echo $CMD
  if [[ ! -f ${OUT_DIR}/${BASE_NAME}.DONE ]]; then
    rm -f ${LOG_DIR}/parse_${BASE_NAME}.stdout
    SGE_CMD="qsub -q long.q -cwd -j y -o ${LOG_DIR}/parse_${BASE_NAME}.stdout -N parse_${BASE_NAME} -S /bin/sh $CMD"
    echo $SGE_CMD
    eval $SGE_CMD
  fi
done
