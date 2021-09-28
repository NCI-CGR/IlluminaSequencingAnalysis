#!/bin/sh

# This script scans and syncronizes the BAM and the BAM_recalibrated folders
# If a BAM file in the BAM folder is newer than the same file in the BAM_recalibrated folder,
# that means a new BAM file is generated and needs recalibration.
# In order to do that, the script spawn a bash script job recalibrate_bam.sh to the cluster


SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")


# source global configuration files
. ${SCRIPT_DIR:-.}/global_config_bash.rc

if [ $# -eq 1 ]; then
     if [[ "$1" == "SOMATIC" ]] || [[ "$1" == "GERMLINE" ]]; then
         VARIANT_TYPE=$1
         MANIFEST=NA
     else
         MANIFEST=$1
         if [[ ! -f $MANIFEST ]]; then
            echo "Error: $MANIFEST is not existed or please specify GERMLINE or SOMATIC as the project type!"
            exit 1
         fi

         VARIANT_TYPE=GERMLINE
     fi
fi
if [ $# -eq 2 ]; then
	# if there is an argument, enter patch mode
	MANIFEST=$1
        if [[ ! -f $MANIFEST ]]; then
            echo "Error: $MANIFEST is not existed!"
            exit 1
        fi 
        VARIANT_TYPE=$2
        if [[ "$VARIANT_TYPE" != "SOMATIC" ]] && [[ "$VARIANT_TYPE" != "GERMLINE" ]]; then
            echo "Error: the second parameter must be SOMATIC or GERMLINE!"
            exit 1
        fi
fi
if [ $# -gt 2 ]; then
        echo "Usage: four forms:"
        echo "1. Only check have any new BAMs in the groups in Manifestsh"
        echo " sh ./step8_sync_and_recalibrate_bam.sh [MANIFEST FILE] GERMLINE|SOMATIC"
        echo " 2. check have any new BAMs in all the BAM_original and processed using the specific variant type"
        echo "sh ./step8_sync_and_recalibrate_bam.sh GERMLINE|SOMATIC"
        echo " 3. check all new BAMs in all the BAM_original with specific project type"
        echo " sh ./step8_sync_and_recalibrate_bam.sh GERMLINE|SOMATIC"
        echo " 4. check have any new BAMs in all the BAM_original and process that as Germline"
        echo " sh ./step8_sync_and_recalibrate_bam.sh"
        echo " If current recalication is for a somatic calling project, will generate additional BAMs with BQSR and without properly paired; If no option specified, will go with GERMLINE"
        exit 1
fi
if [ $# -eq 0 ]; then
       MANIFEST=NA
       VARIANT_TYPE=GERMLINE
fi
if [[ "$MANIFEST" != "NA" ]]; then
	echo "Screening only groups from manifest $MANIFEST_FILE file."
	
	#test if the manifest file is txt file, change to txt if the file is csv
	
	GROUP_COL=$(awk -F "," 'BEGIN{found=0; IGNORECASE = 1} {for (i=1;i<=NF;i++) if ($i=="GROUP") {print i; found=1; exit;}} END {if (found==0) {print 0}}' $MANIFEST)
	DISEASE_GROUPS=""
        if [[ $GROUP_COL -eq 0 ]]; then
           echo "Error: cannot find the GROUP column in the manifest file!"
           exit 1
        fi
	for i in $(awk -F "," -v col=$GROUP_COL '{if (NR>1){print $col}}' $MANIFEST | sort | uniq);do 
	   DISEASE_GROUPS="$DISEASE_GROUPS ${BAM_REFORMATTED_ORIGINAL_DIR}/${i}"
	done
else
	echo "No manifest file. Screening all groups in $BAM_ORIGINAL_DIR folder "
	DISEASE_GROUPS=`ls -d ${BAM_REFORMATTED_ORIGINAL_DIR}/*`
fi


echo "DISEASE_GROUPS: "${DISEASE_GROUPS}
# exit 1
# for GROUP in $BAM_ORIGINAL_DIR/* ; do
for GROUP in $(echo $DISEASE_GROUPS); do
  GROUP_NAME=`basename $GROUP`
  echo $GROUP_NAME
  if [[ -d $GROUP ]]; then  
    for IN_BAM in $GROUP/*.bam; do

      if [[ ! -f $IN_BAM ]]; then
	      continue
      fi

      BAM_NAME=`basename $IN_BAM`
      echo "find ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP_NAME} -name $BAM_NAME -newer $IN_BAM 2>/dev/null | wc -l"
      HAS_NEWER=`find ${BAM_REFORMATTED_RECALIBRATED_DIR}/$GROUP_NAME -name $BAM_NAME -newer $IN_BAM 2>/dev/null | wc -l`
      echo "$BAM_NAME $HAS_NEWER"
      
      # this is a working flag indicating that an existing process
      # is already on this BAM file
      # Move flag files to a up level ${BAM_RECALIBRATED_DIR}  to avoid processing same BAM in different GROUP
      FLG_WORKING=${BAM_REFORMATTED_RECALIBRATED_DIR}/${BAM_NAME}.$RECALIBRATION_WORKING_FLAG_EXTENSION     
      FLG_INQUEUE=${BAM_REFORMATTED_RECALIBRATED_DIR}/${BAM_NAME}.$RECALIBRATION_INQUEUE_FLAG_EXTENSION
      # echo "Flg_working:$FLG_WORKING"
      # echo "Flg_inqueue:$FLG_INQUEUE"

      if [[ $HAS_NEWER -eq 0 ]] && [[  ! -f $FLG_WORKING ]] && [[  ! -f $FLG_INQUEUE ]]; then
        echo
        echo "Launching recalibration job for $BAM_NAME ..."
        if [[ ! -d ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP_NAME} ]]; then
          mkdir -p ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP_NAME}
        fi
        touch $FLG_INQUEUE

	      DATE=`echo $(date +%Y%m%d)`
	      LOG_DIR=${CLUSTER_JOB_LOG_DIR}/${DATE}
        if [[ ! -d $LOG_DIR ]]; then
	        mkdir -p $LOG_DIR
        fi
 	      cd $SCRIPT_DIR
# 	     CMD="qsub -q $QUEUE \
#              -o ${LOG_DIR}/_recalibration_${BAM_NAME}.stdout -e ${LOG_DIR}/_recalibration_${BAM_NAME}.stderr \
#              -N RECAL.$BAM_NAME -v SCRIPT_DIR=$SCRIPT_DIR \
#              -S /bin/sh ${SCRIPT_DIR}/recalibrate_bam.sh $IN_BAM ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP_NAME}/${BAM_NAME} $VARIANT_TYPE"

	      #Job Script in Biowulf -->
        CMD="sbatch --time=10-00:00:00 \
                    --mem=20G \
                    --job-name=RECAL.${BAM_NAME} \
                    --export=SCRIPT_DIR='${SCRIPT_DIR}' \
                    --output=${LOG_DIR}/_recalibration_${BAM_NAME}.stdout \
                    --error=${LOG_DIR}/_recalibration_${BAM_NAME}.stderr \
                    --wrap=\"bash ${SCRIPT_DIR}/recalibrate_bam.sh ${IN_BAM} ${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP_NAME}/${BAM_NAME} ${VARIANT_TYPE}\""
        #<--
	      echo $CMD
	      eval $CMD
      fi
    done
  fi
done
