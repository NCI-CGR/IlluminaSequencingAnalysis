#!/bin/sh

# Copy new incoming BAM files (from BAM_new_incoming) to the original BAM files folder (BAM_original)
# so that they can be recalibrated (into BAM_recalibrated) for future builds
# At the same time, backup old BAM files that are to be replaced by the new files to the backup folder (BAM_old_backups)

# Source global configurations

SCRIPT=$(readlink -f "$0")
SCRIPT_HOME=$(dirname "$SCRIPT")

#. ${SCRIPT_HOME}/global_config_bash.rc
SCRIPT=$(readlink -f "$0")
DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc

#Get Working flag and done flag
strFlagWorking=$1
strFlagDone=$2

echo "strFlagWorking: ${strFlagWorking}"
echo "strFlagDone   : ${strFlagDone}"
echo

# . ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc
# BAM_REFORMATTED_ORIGINAL_DIR
# BAM_REFORMATTED_RECALIBRATED_DIR
module load samtools
for SUBDIR in ${BAM_INCOMING_DIR}/*; do
  SUBDIR_NAME=$(basename $SUBDIR)
  for BAM in ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/*.bam; do
    
    if [[ ! -f $BAM ]]; then
      continue
    fi
    USER=`whoami`  
    COUNT_OWNER=`ls -l $BAM | grep $USER | wc -l`
    if [[ "$COUNT_OWNER" == "0" ]]; then
      continue
    fi
    BAM_NAME=$(basename $BAM .bam)

    BAM_ORIGINAL=${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME}/${BAM_NAME}.bam

    BAM_RECALIBRATED=${BAM_REFORMATTED_RECALIBRATED_DIR}/${SUBDIR_NAME}/${BAM_NAME}.bam
    METRICS_FILE=${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME}/${BAM_NAME}_dedup_metrics.txt
    echo "Incoming: $BAM"
    echo "Original: $BAM_ORIGINAL"
    echo "Recalibrated: $BAM_RECALIBRATED"

    # First we wanted to make sure the destination folder is writable by the script
                # before we make any file operations
    if [[ ! -d ${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME} ]]; then
      mkdir -p ${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME} 
      if [[ $? -ne 0 ]]; then
         echo "Error: mkdir ${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME} is failed!"
         exit 1
      fi
    fi

    DATE=$(date +%Y-%m-%d)
    BAM_ORIGINAL_BAK_DIR=${BAM_BACKUP_DIR}/${DATE}/original
    BAM_RECALIBRATED_BAK_DIR=${BAM_BACKUP_DIR}/${DATE}/recalibrated

    if [[ -f $BAM_RECALIBRATED ]]; then
          echo
          echo "Backing up recalibrated BAM file..."
          CMD="mkdir -p $BAM_RECALIBRATED_BAK_DIR"
          echo $CMD
          eval $CMD
          
          SOFT_LINKs=`find /DCEG/Projects/Exome/builds/build_*/bam_location -mount -lname $BAM_RECALIBRATED`
         
          echo $SOFT_LINKs
          for ss in $SOFT_LINKs; do
             if [[ ! -L $ss ]]; then
                CMD="rm -f $ss ${ss}.bai"
                echo $CMD
                eval $CMD
             fi
          done
          
          SOFT_LINK2s=`find /DCEG/Projects/Exome/builds/build_*/bam_location -mount -lname ${BAM_REFORMATTED_RECALIBRATED_DIR}/${SUBDIR_NAME}/${BAM_NAME}_bqsr_final_out.bam`
          echo $SOFT_LINK2s
          for ss in $SOFT_LINK2s; do
             if [[ ! -L $ss ]]; then
                CMD="rm -f $ss ${ss}.bai"
                echo $CMD
                eval $CMD
             fi
          done                  
          CMD="mv ${BAM_REFORMATTED_RECALIBRATED_DIR}/${SUBDIR_NAME}/${BAM_NAME}*.* $BAM_RECALIBRATED_BAK_DIR"
          echo $CMD
          eval $CMD
        
          if [[ $? -ne 0 ]]; then
            echo "Error: backup recalibrated file is failed!"
            exit 1
          fi

          for ss in $SOFT_LINKs; do
             CMD="ln -s ${BAM_RECALIBRATED_BAK_DIR}/${BAM_NAME}.bam $ss"
             echo $CMD
             eval $CMD
             CMD="ln -s ${BAM_RECALIBRATED_BAK_DIR}/${BAM_NAME}.bam.bai ${ss}.bai"
             echo $CMD
             eval $CMD
          done

          for ss in $SOFT_LINK2s; do
             CMD="ln -s ${BAM_RECALIBRATED_BAK_DIR}/${BAM_NAME}_bqsr_final_out.bam $ss"
             echo $CMD
             eval $CMD
             CMD="ln -s ${BAM_RECALIBRATED_BAK_DIR}/${BAM_NAME}_bqsr_final_out.bam.bai ${ss}.bai"
             echo $CMD
             eval $CMD

          done
  

    fi
    # move old original BAM files to backup folder
    # as new BAM files will soon come in
    if [[ -f $BAM_ORIGINAL ]]; then
      echo
      echo "Backing up original BAM file..."
      CMD="mkdir -p $BAM_ORIGINAL_BAK_DIR"
      echo $CMD
      eval $CMD
      CMD="mv ${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME}/${BAM_NAME}*.* $BAM_ORIGINAL_BAK_DIR"
      echo $CMD
      eval $CMD
  
      if [[ $? -ne 0 ]]; then
        echo "Error: backup original file is failed!"
        exit 1
      fi
    fi

    if [[ -f $METRICS_FILE ]]; then
      echo
      echo "Backing up dedup metrics file ..."
      CMD="mv $METRICS_FILE $BAM_ORIGINAL_BAK_DIR"
      echo $CMD
      eval $CMD
      if [[ $? -ne 0 ]]; then
        echo "Error: backup the dedup metrics file $METRICS_FILE is failed!"
        exit 1
      fi

    fi
    echo
    echo "Moving incoming BAM file to 'original' folder..."
    CMD="mv $BAM ${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME} && chmod g+r $BAM_ORIGINAL"
                 
    echo $CMD
    eval $CMD
    if [[ $? -ne 0 ]]; then
       echo "Error: mv $BAM from incoming to original is failed!"
       exit 1
    fi
    echo
    echo "Moving incoming dedup metric file to original folder ..."
    if [[ -f ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/${BAM_NAME}_metrics.txt ]]; then
    # CMD="mv ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/${BAM_NAME}_dedup_metrics.txt ${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME} && chmod g+r $METRICS_FILE"
      CMD="mv ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/${BAM_NAME}_metrics.txt $METRICS_FILE && chmod g+r $METRICS_FILE"
    else
      CMD="mv ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/${BAM_NAME}_dedup_metrics.txt ${BAM_REFORMATTED_ORIGINAL_DIR}/${SUBDIR_NAME} && chmod g+r $METRICS_FILE"  

    fi
    echo $CMD
    eval $CMD
    if [[ $? -ne 0 ]]; then
       echo "Error: mv ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/${BAM_NAME}_dedup_metrics.txt from incoming to original is failed!"
       exit 1
    fi
                

    # guess the incoming BAM index file name
    # and make sure index are always named *.bam.bai
               
    if [[ -f ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/${BAM_NAME}.bai ]]; then
      echo "Moving incoming BAM index file (BAI) to 'original' folder..."
      CMD="mv ${BAM_INCOMING_DIR}/${SUBDIR_NAME}/${BAM_NAME}.bai ${BAM_ORIGINAL}.bai && chmod g+r ${BAM_ORIGINAL}.bai"
    elif [[ -f ${BAM}.bai ]]; then
      echo "Moving incoming original BAM index file (BAI) to 'original' folder..."
      CMD="mv ${BAM}.bai ${BAM_ORIGINAL}.bai && chmod g+r ${BAM_ORIGINAL}.bai"
    else
      echo "BAM index file (*.bai) not found in incoming folder. Creating one..."
      CMD="samtools index $BAM_ORIGINAL ${BAM_ORIGINAL}.bai && chmod g+r ${BAM_ORIGINAL}.bai"
    fi
    echo $CMD
    eval $CMD
          
    if [[ $? -ne 0 ]]; then
      echo "Error: mv index file is failed!"
      exit 1
    fi

    echo 
  done
done  
echo "All done!"

rm ${strFlagWorking}
touch ${strFlagDone}
