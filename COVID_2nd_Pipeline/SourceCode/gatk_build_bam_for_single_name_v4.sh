#!/bin/sh

# Build the BAM file for analysis ID according to list of sample-level BAM file(s)
# The 1st argument shall be the analysis ID, and the rest argument(s) shall be the full paths of the sample-level BAM files

# Source global configurations
# ./global_config_bash.rc
#. ${SCRIPT_HOME}/global_config_bash.rc
SCRIPT=$(readlink -f "$0")
DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc

set -o pipefail
# set -e
module load samtools/1.13 java/1.8.0_211
ANALYSIS_ID=$1
DOWNSAMPLE_RATIO=$2
MANIFEST=$3
strFlagWorking=$4
strFlagDone=$5
shift
shift
shift
shift
shift
PRE_INPUT_LIST=$*
echo "$(date): Running merging..."
#remove duplicates in input list due to manifest duplication
#INPUT_LIST=$(for i in $(echo $PRE_INPUT_LIST | cut -f1- -d' ' | sort | uniq);do echo $i;done | sort | uniq|cut -d'.' -f1)_HQ_paired_dedup_properly_paired_nophix.bam
INPUT_LIST=$(for i in $(echo $PRE_INPUT_LIST | cut -f1- -d' ' | sort | uniq);do echo $i;done | sort | uniq)
echo "INPUT_LIST: $INPUT_LIST"
SUB_DIR=`echo $ANALYSIS_ID | cut -f1 -d_`
#ANALYSIS_ID=`echo $ANALYSIS_ID | sed -e 's/_U//'`
ANALYSIS_ID=`echo "${ANALYSIS_ID//_NOID}"`
echo "TMP_DIR: ${TMP_DIR}"
OUT_TMP_PREFIX=${TMP_DIR}/${ANALYSIS_ID}.copy_or_merge

OUT_TMP_DEDUP=${OUT_TMP_PREFIX}_dedup.bam
OUT_TMP_DEDUP_IDX=${OUT_TMP_PREFIX}_dedup.bai
OUT_TMP_DEDUP_metrics=${OUT_TMP_PREFIX}_dedup.txt

OUT=${BAM_INCOMING_DIR}/${SUB_DIR}/${ANALYSIS_ID}.bam
OUT_IDX=${BAM_INCOMING_DIR}/${SUB_DIR}/${ANALYSIS_ID}.bai
OUT_METRICS=${BAM_INCOMING_DIR}/${SUB_DIR}/${ANALYSIS_ID}_dedup_metrics.txt

if [[ -f $OUT ]] && [[ -f $OUT_IDX ]] && [[ -f $OUT_METRICS ]]; then
  echo "The merged BAM file of this subject $ANALYSIS_ID exists under ${BAM_INCOMING_DIR}, skip merging."
  exit 1
fi


if [[ ! -d ${BAM_INCOMING_DIR}/${SUB_DIR} ]]; then
   mkdir -p ${BAM_INCOMING_DIR}/${SUB_DIR} 2> /dev/null
fi
#MANIFEST=../Manifest_new_incoming.txt
INPUT_LIST2=""

#do the downsample for certain nextseq data
# if [[ "$DOWNSAMPLE_RATIO" -ne 1 ]]; then

DATE=`date '+%Y-%m-%d-%H-%M'`
for bam in $INPUT_LIST; do
  BAM_NAME_ROOT=$(basename $bam .bam)
  BAM_RAW_DIR=$(dirname $bam)	
  BAM_FULL_NAME=$(ls ${BAM_RAW_DIR}/${BAM_NAME_ROOT}_HQ_paired_dedup*.bam) 
  BAI_FULL_NAME=$(ls ${BAM_RAW_DIR}/${BAM_NAME_ROOT}_HQ_paired_dedup*.bai) 
  if [[ ! -L $bam ]] && [[ ! -e ${BAM_FULL_NAME} ]]; then
    ARCHIVE=$(dirname $bam|cut -d'/' -f4-)
	BAM_ARCHIVE_FILE=$(ls /DCEG_Archive/CGR/${ARCHIVE}/${BAM_NAME_ROOT}_HQ_paired_dedup*.bam)
	BAI_ARCHIVE_FILE=$(ls /DCEG_Archive/CGR/${ARCHIVE}/${BAM_NAME_ROOT}_HQ_paired_dedup*.bai)
    if [[ -e $BAM_ARCHIVE_FILE ]] && [[ -e $BAI_ARCHIVE_FILE ]] ; then
      echo "the lane level BAM file is under DCEG_Archive, not on T-drive, copy such BAM file and its bai file to its original flowcell directory"
      PRE_BAM_DIR=$(ls -d /DCEG_Archive/CGR/${ARCHIVE})
      T_DRIVE=$(echo $PRE_BAM_DIR|cut -d'/' -f4-)
      BAM_DIR=/DCEG/CGF/${T_DRIVE}
	  BAM_TYPE_ROOT=$(echo $BAM_ARCHIVE_FILE|rev|cut -d'/' -f1|cut -d'.' -f2|rev)
	  BAI_TYPE_ROOT=$(echo $BAI_ARCHIVE_FILE|rev|cut -d'/' -f1|cut -d'.' -f2|rev)
	  mkdir -p $BAM_DIR 2> /dev/null
	  cp $BAM_ARCHIVE_FILE ${BAM_DIR}/${BAM_TYPE_ROOT}.bam
	  cp $BAI_ARCHIVE_FILE ${BAM_DIR}/${BAI_TYPE_ROOT}.bai
	    if [[ $? -ne 0 ]]; then
          echo "Error: copy BAM file or bai files from DCEG_Archive to T drive failed!"
          exit 1
        fi
	  ln -s ${BAM_DIR}/${BAM_TYPE_ROOT}.bam ${BAM_DIR}/${BAM_NAME_ROOT}.bam
      ln -s ${BAM_DIR}/${BAI_TYPE_ROOT}.bai ${BAM_DIR}/${BAM_NAME_ROOT}.bai
	elif [[ ! -e $BAM_ARCHIVE_FILE ]]; then
	  echo "Error: the lane-level file $bam is not on T-drive or under DCEG_Archive, please retrieve the files from the retrieve_bam.lst!"
      echo $bam >> /DCEG/Projects/Exome/SequencingData/sync_script_v4/${DATE}_retrieve_bam.lst
	else 
      echo "Error: $BAI_ARCHIVE_FILE does not exist, something went wrong, please do manual check!"
      exit 1
    fi
  elif [[ ! -L $bam ]] && [[ -e ${BAM_FULL_NAME} ]]; then
    echo "The softlink of BAM file is missing on T drive, but the BAM file exists, recreate the softlink of BAM file and its BAI file"
    ln -s $BAM_FULL_NAME $bam
    ln -s $BAI_FULL_NAME ${BAM_RAW_DIR}/${BAM_NAME_ROOT}.bai
  else
    echo "the softlink of the $bam exists, proceed further."
  fi
done

for bam in $INPUT_LIST; do
        echo
        echo yes $bam
        # For this BAM, check if the bam in the old version or not; then replace all RG
        CC=`samtools view -H $bam | awk -F"\t" '$1~/@RG/ {for (i=1;i<=NF;i++) { if ($i=="LB:N/A") {print $i}}}' | wc -l`  
        if [[ $CC -eq 1 ]]; then
        # Okd version of lane level BAM 
           BAM_DIR=`dirname $bam`
           BASE_NAME=`basename $bam .bam`
           OUT_BAM=${BAM_DIR}/${BASE_NAME}_reformated.bam
           OLD_NMs=`samtools view -H $bam | awk -F"\t" '$1~/@RG/ {for (i=1;i<=NF;i++) { if ($i~/^ID/) split($i,aa,":"); if ($i~/^PU/) split($i,bb,":"); }} END {print aa[2]":"bb[2]}'`
           OLD_ID=`echo $OLD_NMs | cut -d':' -f1`
           OLD_PU=`echo $OLD_NMs | cut -d':' -f2`
           FLOWCELL_ID=`echo $bam | rev | cut -d"/" -f 3 | rev | cut -d"_" -f4`
           LANE=`echo $bam | rev | cut -d"/" -f 1 | cut -d"." -f2 | head -c 1`
           NEW_ID=$FLOWCELL_ID"."${LANE}
           SED_CMD="s/\bRG:Z:${OLD_ID}\t\b/RG:Z:${NEW_ID}\t/g; s/\b\tPU:Z:${OLD_PU}\b//g; s/\tLB:Z:N\/A//g;"
           CMD="samtools view -h $bam | sed '$SED_CMD' | samtools view -S -b - > ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam"

           echo $CMD
           eval $CMD
           if [[ $? -ne 0 ]]; then
             echo "Error: samtools is failed to modify the body part!"
             exit 1
           fi

           SEQ_DATE=`echo $bam | rev | cut -d"/" -f 3 | rev | cut -d"_" -f1`
           INDEX=`echo $bam | rev | cut -d"/" -f 1 | cut -d"_" -f2 | rev`
           # SAMPLE_ID=`echo $bam | rev | cut -d"/" -f 1 | rev | cut -d"_" -f1`
           SAMPLE_ID=`echo $bam |rev | cut -d"/" -f 1 |cut -d'_' -f3-|rev`
           NEW_PU=$FLOWCELL_ID"20"$SEQ_DATE"."$LANE"."$INDEX
           NEW_LB=${SAMPLE_ID}"_"${INDEX}
           
           CMD="samtools view -H ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam | awk -F\"\t\" -v id=\"$NEW_ID\" -v lb=\"$NEW_LB\" -v pu=\"$NEW_PU\" '{if (\$1~/^@RG/){ 
                                                for (i=1;i<=NF;i++) {
                                                     if (\$i~/SM:/){
                                                       sm=\$i;
                                                       break;
                                                     }
                                                  }
                                                  print \"@RG\tID:\"id\"\t\"sm\"\tLB:\"lb\"\tPL:ILLUMINA\tPU:\"pu\"\tCN:CGR\";
                                                }
                                                else
                                                  print
                                                 }' > ${TMP_DIR}/${ANALYSIS_ID}_header.sam"
           echo
           echo $CMD
           eval $CMD
           if [[ $? -ne 0 ]]; then
             echo "Error: samtools is failed to modify the header part!"
             exit 1
           fi

           CMD="samtools reheader ${TMP_DIR}/${ANALYSIS_ID}_header.sam ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam >  ${TMP_DIR}/${ANALYSIS_ID}_tmp2.bam"
           echo
           echo $CMD
           eval $CMD

           if [[ $? -ne 0 ]]; then
              echo "Error: samtools reheader is failed!"
              exit 1
           fi
           echo
           CMD="samtools sort -@ 8 -T ${TMP_DIR}/${ANALYSIS_ID}_temp -o ${TMP_DIR}/${ANALYSIS_ID}_sorted.bam ${TMP_DIR}/${ANALYSIS_ID}_tmp2.bam"
           echo $CMD
           eval $CMD
           if [[ $? -ne 0 ]]; then
              echo "Error: samtools sorting is failed!"
              exit 1
           fi
           echo "$(date) samtools sorting bam was finished successfully!"
           #CMD="mv ${TMP_DIR}/${ANALYSIS_ID}_sorted.bam $OUT_BAM && samtools index $OUT_BAM ${BAM_DIR}/${BASE_NAME}_reformated.bai"
		   CMD="cp ${TMP_DIR}/${ANALYSIS_ID}_sorted.bam $OUT_BAM && samtools index $OUT_BAM ${BAM_DIR}/${BASE_NAME}_reformated.bai"
           echo
           echo $CMD
           eval $CMD
           if [[ $? -ne 0 ]]; then
              echo "Error: $(date) sam indexing was failed!"
              exit 1
           else
             echo "$(date) sam indexing was finished successfully!"
           fi
           OLD_COUNT=`samtools view -@ 8 -c $bam`
           NEW_COUNT=`samtools view -@ 8 -c $OUT_BAM`
           echo $OLD_COUNT" "$NEW_COUNT
           if [[ $OLD_COUNT -ne $NEW_COUNT ]]; then
             echo "Error: the count is inconsistent!"
             echo
             echo "Cleaning up ..."
             CMD="rm -f $OUT_BAM ${TMP_DIR}/${ANALYSIS_ID}_header.sam ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam ${TMP_DIR}/${ANALYSIS_ID}_tmp2.bam ${TMP_DIR}/${ANALYSIS_ID}_temp* ${TMP_DIR}/${ANALYSIS_ID}_sorted.bam"
             echo $CMD
             eval $CMD
             exit 1
           fi
           echo 
           echo "Cleaning up ..." 
           CMD="rm -f $bam ${BAM_DIR}/${BASE_NAME}.bai ${TMP_DIR}/${ANALYSIS_ID}_header.sam ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam ${TMP_DIR}/${ANALYSIS_ID}_tmp2.bam ${TMP_DIR}/${ANALYSIS_ID}_temp* ${TMP_DIR}/${ANALYSIS_ID}_sorted.bam"
           echo $CMD
           eval $CMD
           if [[ $? -ne 0 ]]; then
              echo "Error: remove softlink is failed!"
              exit 1
           fi 
           CMD="ln -s $OUT_BAM $bam"
           echo $CMD
           eval $CMD
           if [[ $? -ne 0 ]]; then
              echo "Error: cannot create softlink!"
              exit 1
           fi
           CMD="ln -s ${BAM_DIR}/${BASE_NAME}_reformated.bai ${BAM_DIR}/${BASE_NAME}.bai"
           echo $CMD
           eval $CMD
           if [[ $? -ne 0 ]]; then
              echo "Error: cannot create softlink!"
              exit 1
           fi
        else
          # New bam or LB:N/A greater than 1
          if [[ $CC -gt 1 ]]; then
            echo "Error: there are more than one LB in the @RG of $bam"
            exit 1
          fi
        fi

        if [[ $bam == *"NextSeq"*  ]]; then
           if [[ "$DOWNSAMPLE_RATIO" != "1" ]]; then
               echo "[$(date)] Downsample for $bam"
               CMD="java -Xmx8g -jar $PICARD DownsampleSam \
                        --INPUT $bam \
                        --OUTPUT ${TMP_DIR}/${ANALYSIS_ID}.downsample.bam \
                        --PROBABILITY $DOWNSAMPLE_RATIO"
               echo $CMD
               eval $CMD
               if [[ $? -ne 0 ]]; then
                   echo "Error: $(date) Downsampling is failed!"
                   exit 1
               else
                   echo "[$(date)] Downsampling is done!"
               fi
               INPUT_LIST2="$INPUT_LIST2 ${TMP_DIR}/${ANALYSIS_ID}.downsample.bam"
           else
                INPUT_LIST2="$INPUT_LIST2 $bam"
           fi
        else
           INPUT_LIST2="$INPUT_LIST2 $bam"
        fi
done
# fi

#echo $DOWNSAMPLE_RATIO
echo ====================
date
echo Building BAM file for $ANALYSIS_ID
echo Input: $INPUT_LIST
echo Downsampled Input: $INPUT_LIST2
echo Output: $OUT
echo ====================

if [[ -f $OUT ]]; then
   echo "Error: The output file $OUT already exists. Operation skipped."
   exit 1
fi
echo "Multiple input BAMs. cleaning up..."
CMD="rm -f ${OUT_TMP_PREFIX}1*.* ${OUT_TMP_PREFIX}.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
     echo "Error: rm temp files failed!"
     exit 1
fi


if [[ $# -eq 0 ]]; then
   echo "Invalid argument. Expecting at least one input BAM file." 1>&2
   exit 1
elif [[ $# -eq 1 ]]; then
   #CMD="cp -l $INPUT_LIST ${OUT_TMP_PREFIX}1.bam"
   #CMD="mv $INPUT_LIST ${OUT_TMP_PREFIX}1.bam"
   CMD="cp $INPUT_LIST ${OUT_TMP_PREFIX}1.bam"
   
   echo $CMD
   #eval $CMD
   cp ${INPUT_LIST} ${OUT_TMP_PREFIX}1.bam

   if [[ $? -ne 0 ]]; then
     echo "Copy Failed!"
   else
     echo "Copy Successfully!"
   fi
else
   echo "[$(date)]: samtools merging ... "
   CMD="samtools merge -@ 8 -f ${OUT_TMP_PREFIX}1.bam $INPUT_LIST2"
   echo $CMD
   eval $CMD
   if [[ $? -ne 0 ]]; then
       echo "Error: samtools merge is failed!"
       exit 1
   fi
   echo "[$(date)]: Merging is done!"
fi

echo
echo "Copy INPUT_LIST finished!"

MERGED_COUNT=`samtools view -@ 8 -c ${OUT_TMP_PREFIX}1.bam`

echo
echo "[$(date)]: Starting reordersam using Picard ..."
CMD="java -Xmx16g -jar $PICARD ReorderSam \
          --INPUT ${OUT_TMP_PREFIX}1.bam \
          --OUTPUT ${OUT_TMP_PREFIX}.bam \
          --SEQUENCE_DICTIONARY ${REFERENCE_GENOME} \
          --REFERENCE_SEQUENCE ${REFERENCE_GENOME} \
          --VALIDATION_STRINGENCY STRICT \
          --TMP_DIR $TMP_DIR"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
     echo "Error: Picard reordersam is failed!"
     exit 1
fi

if [[ ! -f ${OUT_TMP_PREFIX}.bam ]]; then
    echo "Error: the expected output file ${OUT_TMP_PREFIX}.bam is missing!"
    exit 1
fi
echo
# echo "[$(date)]: Removing SM:i:* in the body of ${OUT_TMP_PREFIX}.bam ..."
# CMD="samtools view -h ${OUT_TMP_PREFIX}.bam | sed 's/\tSM:i:*//g' | samtools view -S -b - > ${OUT_TMP_PREFIX}_SM_removed.bam"
# echo $CMD
# eval $CMD

# if [[ $? -ne 0 ]]; then
#    echo "Error: SM:i is failed to be removed!"
#    exit 1
# fi

# if [[ ! -f ${OUT_TMP_PREFIX}_SM_removed.bam ]]; then
#    echo "Error: expected output file ${OUT_TMP_PREFIX}_SM_removed.bam does not exist!"
#    exit 1
# fi
echo "[$(date)]: Picard reordersam ${OUT_TMP_PREFIX}.bam is done!"

echo
echo "[$(date)]: Reheader SM to the sample name in the @RG part ...."
CMD="samtools view -H ${OUT_TMP_PREFIX}.bam | awk -F\"\t\" -v sm=\"$ANALYSIS_ID\" -v OFS=\"\t\" '{if (\$1~/^@RG/) {
                                                for (i=1;i<=NF;i++){
                                                  if (\$i~/SM:/)
                                                    \$i=\"SM:\"sm;
                                                }
                                                print;
                                                }
                                                else
                                                  print;
                                                }' > ${OUT_TMP_PREFIX}_SM_renamed.sam"
echo
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: samtools is failed to output the header part!"
   exit 1
fi

CMD="samtools reheader ${OUT_TMP_PREFIX}_SM_renamed.sam ${OUT_TMP_PREFIX}.bam > ${OUT_TMP_PREFIX}_SM_renamed.bam"
echo
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: samtools is failed to reheader part!"
   exit 1
fi
if [[ ! -f ${OUT_TMP_PREFIX}_SM_renamed.bam ]]; then
    echo "Error: the output BAM: ${OUT_TMP_PREFIX}_SM_renamed.bam is missing!"
    exit 
fi
echo "[$(date)]: Reheader SM is done!"
echo
echo "[$(date)]: Sorting the merged file ${OUT_TMP_PREFIX}_SM_renamed.bam ..."
CMD="samtools sort -@ 8 -T ${OUT_TMP_PREFIX}_SM_renamed_sorted_temp -o ${OUT_TMP_PREFIX}_SM_renamed_sorted.bam ${OUT_TMP_PREFIX}_SM_renamed.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
    echo "Error: samtools sorting is failed!"
    exit 1
fi


echo
echo "The first checkpoint: check if the read counts are consistent"
SECOND_MERGED_COUNT=`samtools view -@ 8 -c ${OUT_TMP_PREFIX}_SM_renamed_sorted.bam`
if [[ $SECOND_MERGED_COUNT -ne $MERGED_COUNT ]]; then
   echo "Error: the counts: $SECOND_MERGED_COUNT and $MERGED_COUNT are inconsistent!"
   exit 1
else
   echo "Read count: $MERGED_COUNT"
fi

echo
echo "The second checkpoint: check if the PU and LB are consistent with MANIFEST (in case the error in flowcell folernames or file names)"
A_ID=`echo $ANALYSIS_ID | cut -d_ -f2-`
echo "A_ID: "${A_ID}
RG=`awk -F"," -v analysis_id=$A_ID 'BEGIN{i=0}
 {
  if ($13==analysis_id){
    SEQ_DATE=$2
    split($2,aa,"/")
    if (length(aa[1])==1)
      aa[1]=0aa[1]
    if (length(aa[2])==1)
      aa[2]=0aa[2]
    RG_PUs[i]=$3aa[3]aa[1]aa[2]"."$4"."$5
    RG_LBs[i]=$6"_"$5
    i++
  }
}END{
   for (i=0;i<length(RG_LBs);i++){
     if (i < length(RG_LBs)-1)
       printf RG_LBs[i]"|"RG_PUs[i]";"
     else
       printf RG_LBs[i]"|"RG_PUs[i]
   }

} ' $MANIFEST`


BAM_RG=`samtools view -H ${OUT_TMP_PREFIX}_SM_renamed_sorted.bam | awk -F"\t" 'BEGIN{j=0}
{
   if ($1~/^@RG/){
     for (i=1;i<=NF;i++){
        if ($i~/^LB:/){
          split($i,aa,":");
          RG_LBs[j]=aa[2];
        }
        if ($i~/^PU:/) {
          split($i,aa,":");
          RG_PUs[j]=aa[2];
        }
     }
#     if (RG_IDs[j]!=RG_PUs[j]){
#       printf "Error: ID: "RG_IDs[j]" is not equal to PU:" RG_PUs[j]
#       exit 1
#     }
     j++
   }
}END{
   for (i=0;i<length(RG_LBs);i++){
     if (i < length(RG_LBs)-1)
       printf RG_LBs[i]"|"RG_PUs[i]";"
     else
       printf RG_LBs[i]"|"RG_PUs[i]
   }
}'`

IFS='; ' read -r -a RG_M <<< "$RG"
IFS='; ' read -r -a RG_B <<< "$BAM_RG"
M_LEN=${#RG_M[@]}
B_LEN=${#RG_B[@]}
echo
echo "M_LEN: "${M_LEN}
echo "B_LEN: "${B_LEN}
echo
if [[ $B_LEN -ne $M_LEN ]]; then
   HISEQ_NUM=`awk -F"," -v analysis_id=$A_ID 'BEGIN{dd=0; hiseq_num=0; i=0}
   {
    if ($13==analysis_id){
      if($1!~/^E0296/){ 
        dd+=$11;
        hiseq_num++;
      }
    }
   }END{
      if (dd>=60)
        printf hiseq_num
      else
        printf dd
        printf "Err"
   }' $MANIFEST` 
   echo "HiSeq Lane Numbers: $HISEQ_NUM"
   if [[ $HISEQ_NUM == "Err" ]]; then
     echo "Error: the total depth of hiseq lanes in the $MANIFEST $M_LEN less than 60!"
     exit 1
   else
     if [[ "$HISEQ_NUM" != "$B_LEN" ]]; then
       echo "Error: the number of lanes in the BAM ${OUT_TMP_PREFIX}_SM_renamed_sorted.bam $B_LEN is diffrent than the number in the $MANIFEST $M_LEN !"
       exit 1
     else
       # Now can skip nextseq lane 
       RG=`awk -F"," -v analysis_id=$A_ID 'BEGIN{i=0}
       {
        if (($13==analysis_id) && ($1!~/^E0296/)) {
            SEQ_DATE=$2
            split($2,aa,"/")
            if (length(aa[1])==1)
              aa[1]=0aa[1]
            if (length(aa[2])==1)
               aa[2]=0aa[2]
            RG_PUs[i]=$3aa[3]aa[1]aa[2]"."$4"."$5
            RG_LBs[i]=$6"_"$5
            i++
        }
       }END{
         for (i=0;i<length(RG_LBs);i++){
           if (i < length(RG_LBs)-1)
             printf RG_LBs[i]"|"RG_PUs[i]";"
           else
             printf RG_LBs[i]"|"RG_PUs[i]
         }
        }' $MANIFEST`
        IFS='; ' read -r -a RG_M <<< "$RG"
     fi
   fi
fi


for i in "${!RG_M[@]}"; do
   FOUND=0
   for j in "${!RG_B[@]}"; do
     echo ${RG_M[i]}" vs "${RG_B[j]}
     if [[ "${RG_M[i]}" == "${RG_B[j]}" ]]; then
       FOUND=1
       break
     fi
   done
   if [[ $FOUND -eq 0 ]]; then
      echo "Error: ${RG_M[i]} cannot find the matches in ${OUT_TMP_PREFIX}_SM_renamed_sorted.bam!"
      exit 1
   fi
done

echo

echo "[$(date)]: Remove duplicates for '$ANALYSIS_ID'..."
CMD="java -Xmx8g -jar $PICARD MarkDuplicates \
          --INPUT ${OUT_TMP_PREFIX}_SM_renamed_sorted.bam \
          --OUTPUT $OUT_TMP_DEDUP \
          --METRICS_FILE $OUT_TMP_DEDUP_metrics
          --ASSUME_SORTED true \
          --REMOVE_DUPLICATES true \
          --VALIDATION_STRINGENCY STRICT \
          --CREATE_INDEX true"
echo $CMD
eval $CMD

if [[ $? -ne 0 ]]; then
    echo "Error: $(date) dedup is failed!"
    exit 1
fi


if [[ ! -f ${OUT_TMP_DEDUP} ]]; then
    echo "Error: expected output file $OUT_TMP_DEDUP does not exist!"
    exit 1
fi

echo "[$(date)]: Doing the last piece of housekeeping..."

#CMD="mv $OUT_TMP_DEDUP $OUT && mv $OUT_TMP_DEDUP_IDX $OUT_IDX && mv $OUT_TMP_DEDUP_metrics $OUT_METRICS"
CMD="cp $OUT_TMP_DEDUP $OUT && cp $OUT_TMP_DEDUP_IDX $OUT_IDX && cp $OUT_TMP_DEDUP_metrics $OUT_METRICS"
echo $CMD
eval $CMD
echo
CMD="rm -f ${OUT_TMP_PREFIX}.bam ${OUT_TMP_PREFIX}.bam.bai $OUT_TMP_DEDUP ${OUT_TMP_DEDUP}.bai $OUT_TMP_DEDUP_IDX ${OUT_TMP_PREFIX}.bai ${OUT_TMP_PREFIX}1.bam $OUT_TMP_DEDUP_metrics ${TMP_DIR}/${ANALYSIS_ID}.downsample.bam ${OUT_TMP_PREFIX}_SM_removed.bam  ${OUT_TMP_PREFIX}_SM_renamed.sam ${OUT_TMP_PREFIX}_SM_renamed.bam ${OUT_TMP_PREFIX}_SM_renamed_sorted.bam ${OUT_TMP_PREFIX}_SM_renamed_sorted_temp*"
echo $CMD
eval $CMD

echo "[$(date)] Done!"
rm ${strFlagWorking}
touch ${strFlagDone}
