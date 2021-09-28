#!/bin/sh
# Downsampling the lane-level BAMs
BAM=$1
DOWNSAMPLE_RATIO=$2
SAMPLE=$3
touch ${LOG_DIR}/${SAMPLE}.WORKING

PICARD=/DCEG/Resources/Tools/Picard/Picard-2.18.11/picard.jar
LANE_DIR=`dirname $BAM`
echo $BAM
BASE_NAME=`basename $BAM _HQ_paired_dedup_properly_paired_nophix.bam`
echo $BASE_NAME
CMD="rm -f ${LANE_DIR}/${BASE_NAME}_downsampled.bam ${LANE_DIR}/${BASE_NAME}_temp*.bam ${LANE_DIR}/${BASE_NAME}_downsampled_sorted.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: remove intermeidate files is failed!"
   exit 1
fi
if [[ -f ${LANE_DIR}/${BASE_NAME}.bam ]]; then
  echo "Error: ${LANE_DIR}/${BASE_NAME}.bam is not a softlink!"
  exit 1
fi 

if [[ -L ${LANE_DIR}/${BASE_NAME}.bai ]]; then
  echo "Error: ${LANE_DIR}/${BASE_NAME}.bai is not a softlink!"
  exit 1
fi

if [[ -L ${LANE_DIR}/${BASE_NAME}.bam ]]; then
  rm -f ${LANE_DIR}/${BASE_NAME}.bam
  if [[ $? -ne 0 ]]; then
   echo "Error: remove softlink is failed!"
   exit 1
  fi
fi

if [[ -L ${LANE_DIR}/${BASE_NAME}.bai ]]; then
  rm -f ${LANE_DIR}/${BASE_NAME}.bai
  if [[ $? -ne 0 ]]; then
   echo "Error: remove softlink is failed!"
   exit 1
  fi
fi
 
echo "[$(date)] Downsample for $BAM"
CMD="java -Xmx8g -jar $PICARD DownsampleSam \
                     INPUT=$BAM          \
                     OUTPUT=${LANE_DIR}/${BASE_NAME}_downsampled.bam \
                     PROBABILITY=$DOWNSAMPLE_RATIO"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: $(date) Downsampling is failed!"
   exit 1
else
   echo "[$(date)] Downsampling is done!"
fi
module load samtools
CMD="samtools sort -T ${LANE_DIR}/${BASE_NAME}_temp -o ${LANE_DIR}/${BASE_NAME}_downsampled_sorted.bam ${LANE_DIR}/${BASE_NAME}_downsampled.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: $(date) samtools sorting is failed!"
   exit 1

fi
CMD="rm -f ${LANE_DIR}/${BASE_NAME}_downsampled.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: remove symblick link is failed!"
   exit 1
fi 
CMD="ln -s ${LANE_DIR}/${BASE_NAME}_downsampled_sorted.bam ${LANE_DIR}/${BASE_NAME}.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: create symblick link is failed!"
   exit 1
fi

rm ${LOG_DIR}/${SAMPLE}.WORKING
touch ${LOG_DIR}/${SAMPLE}.DONE
echo "$(date): All done!"
