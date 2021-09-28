#!/bin/sh


SCRIPT=$(readlink -f "$0")
SCRIPT_HOME=$(dirname "$SCRIPT")

. ${SCRIPT_HOME}/global_config_bash.rc
# . ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc
if [[ $# -ne 2 ]]; then
   echo "Usage: sh ./step7a_reformat_bams.sh MANIFEST RESTORE_FILE"
   exit 1
fi
MANIFEST=$1
RESTORE=$2

if [[ -f $RESTORE ]]; then
  rm -f $RESTORE
fi
LOG_DIR=${BAM_REFORMATTED_ORIGINAL_DIR}/LOGs
SAMPLEs=`awk -F"," 'NR>1 {print $7":"$13}' $MANIFEST`
# echo $SAMPLEs
FOUND=0
RES=0
for SAMPLE in $SAMPLEs; do
  # echo $SAMPLE
  GROUP=`echo $SAMPLE | cut -d":" -f1`
  ANALYSIS_ID=`echo $SAMPLE | cut -d":" -f2`
  BAM_NAME="${GROUP}_${ANALYSIS_ID}"
  # echo ${BAM_REFORMATTED_ORIGINAL_DIR}/${GROUP}/${BAM_NAME}.bam 
  # echo ${BAM_INCOMING_DIR}/${GROUP}/${BAM_NAME}.bam 
  # echo ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${BAM_NAME}.bam
  if [[ ! -f ${BAM_REFORMATTED_ORIGINAL_DIR}/${GROUP}/${BAM_NAME}.bam && ! -f ${BAM_INCOMING_DIR}/${GROUP}/${BAM_NAME}.bam && ! -f ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${BAM_NAME}.bam && ! -f ${LOG_DIR}/fix_${ANALYSIS_ID}.WORKING && ! -f ${LOG_DIR}/fix_${ANALYSIS_ID}.INQUEUE && ! -f ${LOG_DIR}/fix_${ANALYSIS_ID}.DONE ]]; then
    if [[ -f ${BAM_ORIGINAL_DIR}/${GROUP}/${BAM_NAME}.bam ]]; then
      FOUND=1     
      rm -f ${LOG_DIR}/fix_${ANALYSIS_ID}.std???
      CMD="qsub -q seq*.q,long.q -o ${LOG_DIR}/fix_${ANALYSIS_ID}.stdout -e ${LOG_DIR}/fix_${ANALYSIS_ID}.stderr -N ${ANALYSIS_ID}_fix -v SCRIPT_HOME=$SCRIPT_HOME,LOG_DIR=$LOG_DIR -S /bin/sh ${SCRIPT_HOME}/misc_fix_bam_header_dedup.sh ${BAM_ORIGINAL_DIR}/${GROUP}/${BAM_NAME}.bam ${BAM_INCOMING_DIR}/${GROUP}/${BAM_NAME}.bam $GROUP $ANALYSIS_ID $MANIFEST ${BAM_INCOMING_DIR}/${GROUP}/${BAM_NAME}_dedup_metrics.txt"
      echo $CMD
      eval $CMD
      touch ${LOG_DIR}/fix_${ANALYSIS_ID}.INQUEUE
      # if [[ "$BAM_NAME" == C* ]]; then
      #  exit 1
      # fi
    else
      RES=1
      echo ${BAM_ORIGINAL_DIR}/${GROUP}/${BAM_NAME}.bam >> $RESTORE
      echo $CMD       
    fi
 # else
 #   if [[ -f ${BAM_REFORMATTED_ORIGINAL_DIR}/${GROUP}/${BAM_NAME}.bam ]]; then
 #     echo "Found ${BAM_REFORMATTED_ORIGINAL_DIR}/${GROUP}/${BAM_NAME}.bam!"
 #   fi
 #   if [[ -f ${BAM_INCOMING_DIR}/${GROUP}/${BAM_NAME}.bam ]]; then
 #     echo "Found ${BAM_INCOMING_DIR}/${GROUP}/${BAM_NAME}.bam!"
 #   fi
 #   if [[ -f ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${BAM_NAME}.bam ]]; then
 #     echo "Found ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${BAM_NAME}.bam!"
 #   fi 
  fi
done
if [[ $FOUND -eq 0 && $RES -eq 0 ]]; then
  echo "All the BAMs are new!"
else
  if [[ $RES -eq 1 ]]; then
    echo "Please check $RESTORE to restore BAMs!"
  fi
fi


