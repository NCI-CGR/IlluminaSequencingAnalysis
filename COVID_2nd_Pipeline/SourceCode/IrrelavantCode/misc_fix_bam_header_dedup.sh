#!/bin/sh
. ${SCRIPT_HOME}/global_config_bash.rc
IN_BAM=$1
OUT_BAM=$2
GROUP=$3
ANALYSIS_ID=$4
MANIFEST=$5
OUT_DEDUP_metrics=$6
# set -e
set -o pipefail
# TMP_DIR=/ttemp/wmy
rm -f ${LOG_DIR}/fix_${ANALYSIS_ID}.INQUEUE
touch ${LOG_DIR}/fix_${ANALYSIS_ID}.WORKING
# IN_BAM=~/tmp2/tmp.bam
# OUT_BAM=~/tmp2/TC/TC_TC_5360_1001_A.bam
# GROUP=TC
# ANALYSIS_ID=TC_5360_1001_A
# MANIFEST=/DCEG/Projects/Exome/builds/2016-09-07_with_population_control/manifest/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv

# IN_BAM=~/tmp2/tmp_WM.bam
# OUT_BAM=~/tmp2/LPD/LPD_WM_4954_1003_A.bam
# GROUP=LPD
# ANALYSIS_ID=WM_4954_1003_A

echo "Cleaning up the tmp folder..."
CMD="rm -f ${TMP_DIR}/${ANALYSIS_ID}*"
echo $CMD
eval $CMD
CMD="rm -f $OUT_BAM"
echo $CMD
eval $CMD

DIR=`dirname $OUT_BAM`
if [[ ! -d $DIR ]]; then
  mkdir -p $DIR
  if [[ $? -ne 0 ]]; then
    echo "Error: failed to create $DIR !"
    exit 1
  fi
fi
module load samtools/1.8
START=`date +%s`
echo "[$(date)]: Starting processing ..."

RG=`awk -F"," -v analysis_id=$ANALYSIS_ID 'BEGIN{i=0}
 {
  if ($13==analysis_id){
    RG_IDs[i]=$3"."$4
    SEQ_DATE=$2
    RG_OLDIDs[i]=$6"_"$5"_L00"$4"_HQ_paired"
    split($2,aa,"/")
    if (length(aa[1])==1)
      aa[1]=0aa[1]
    if (length(aa[2])==1)
      aa[2]=0aa[2]    
    RG_PUs[i]=$3aa[3]aa[1]aa[2]"."$4"."$5
    RG_LBs[i]=$6"_"$5
    RG_INSTRUMENTs[i]=$1
    i++   
  }
}END{
   for (i=0;i<length(RG_IDs);i++){
     if (i < length(RG_IDs)-1)
       printf RG_IDs[i]"|"RG_PUs[i]"|"RG_LBs[i]"|"RG_OLDIDs[i]"|"RG_INSTRUMENTs[i]";"
     else
       printf RG_IDs[i]"|"RG_PUs[i]"|"RG_LBs[i]"|"RG_OLDIDs[i]"|"RG_INSTRUMENTs[i]
   }

} ' $MANIFEST`
if [[ $? -ne 0 ]]; then
   echo "Error: retriving RG from $MANIFEST is failed!"
   exit 1
fi

echo "RG got from $MANIFEST: $RG"
echo
SM=${GROUP}_${ANALYSIS_ID}

OLD_RG=`samtools view -H $IN_BAM | awk -F"\t" 'BEGIN{j=0}
{
   if ($1~/^@RG/){
     for (i=1;i<=NF;i++){
        if ($i~/^ID:/){
          split($i,aa,":");
          RG_IDs[j]=aa[2];
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
   for (i=0;i<length(RG_IDs);i++){
     if (i < length(RG_IDs)-1)
       printf RG_IDs[i]"|"RG_PUs[i]";"
     else
       printf RG_IDs[i]"|"RG_PUs[i]
   }
}'`

if [[ $? -ne 0 ]]; then
   echo "Error: retriving RG from $IN_BAM is failed!"
   exit 1
fi

echo "RG got from $IN_BAM: $OLD_RG"
echo
if [[ $OLD_RG == "Error"* ]]; then
  echo $OLD_RG
  exit 1
fi
IFS='; ' read -r -a OLD <<< "$OLD_RG"
IFS='; ' read -r -a NEW <<< "$RG"
OLD_LEN=${#OLD[@]}
NEW_LEN=${#NEW[@]}
echo "The number of lanes in the old BAM $IN_BAM is: $OLD_LEN and the number in the $MANIFEST is $NEW_LEN"
if [[ $OLD_LEN -ne $NEW_LEN ]]; then
   NEXT_COUNT=0
   for i in "${!NEW[@]}"; do
      IFS='| ' read -r -a NEW_RGs <<< "${NEW[i]}"  
      if [[ "${NEW[4]}" == "E0296"* ]]; then
         NEXT_COUNT=$((NEXT_COUNT+1))
      fi
   done
   if [[ $NEXT_COUNT -eq 1 ]]; then
      echo "Warning: there is one lane in $MANIFEST is nextseq. The number of lanes in the old BAM $IN_BAM $OLD_LEN is diffrent than the number in the $MANIFEST $NEW_LEN!"
   else
      echo "Error: the number of lanes in the old BAM $IN_BAM $OLD_LEN is diffrent than the number in the $MANIFEST $NEW_LEN !"
      exit 1
   fi
fi

SED_CMD1=""
SED_CMD2=""
for i in "${!OLD[@]}"; do
   FOUND=0
   IFS='| ' read -r -a OLD_RGs <<< "${OLD[i]}"
   for j in "${!NEW[@]}"; do
     IFS='| ' read -r -a NEW_RGs <<< "${NEW[j]}"
     echo ${OLD_RGs[0]}" vs "${NEW_RGs[3]}
     if [[ "${OLD_RGs[0]}" == "${NEW_RGs[3]}" ]]; then
       ACTUAL_ID=`samtools view $IN_BAM | awk -F"\t" -v id=${OLD_RGs[0]} '$15=="RG:Z:"id {split($1,aa,":"); print aa[3]"."aa[4]; exit;}'`
       NEW_RGID=`echo ${NEW_RGs[0]} | awk '{print substr($1,2);}'`
       if [[ "$ACTUAL_ID" == "${NEW_RGs[0]}" || "$ACTUAL_ID" == "$NEW_RGID" ]]; then
         FOUND=1
         if [[ ${OLD_RGs[1]} == "N/A" ]]; then
           SED_CMD2="${SED_CMD2}s/\bRG:Z:${OLD_RGs[0]}\t\b/RG:Z:${NEW_RGs[0]}\t/g; s/\b\tPU:Z:[^ ]*//g; "
         else
           SED_CMD2="${SED_CMD2}s/\bRG:Z:${OLD_RGs[0]}\t\b/RG:Z:${NEW_RGs[0]}\t/g; s/\b\tPU:Z:${OLD_RGs[1]}\b//g; "
         fi
         break
       fi
     fi
   done
   if [[ $FOUND -eq 0 ]]; then
      ACTUAL_ID=`samtools view $IN_BAM | awk -F"\t" -v id=${OLD_RGs[0]} '$15=="RG:Z:"id {split($1,aa,":"); print aa[3]"."aa[4]; exit;}'`
      for j in "${!NEW[@]}"; do
        IFS='| ' read -r -a NEW_RGs <<< "${NEW[j]}"
        NEW_RGID=`echo ${NEW_RGs[0]} | awk '{print substr($1,2);}'`
        echo ${OLD_RGs[0]}" "$ACTUAL_ID" vs "${NEW_RGs[0]}" or "$NEW_RGID
        if [[ "$ACTUAL_ID" == "${NEW_RGs[0]}" || "$ACTUAL_ID" == "$NEW_RGID" ]]; then
          FOUND=1
          if [[ ${OLD_RGs[1]} == "N/A" ]]; then
             SED_CMD1="${SED_CMD1}s/\bRG:Z:${OLD_RGs[0]}\t\b/RG:Z:${NEW_RGs[0]}\t/g; s/\b\tPU:Z:[^ ]*//g; "
          else
             SED_CMD1="${SED_CMD1}s/\bRG:Z:${OLD_RGs[0]}\t\b/RG:Z:${NEW_RGs[0]}\t/g; s/\b\tPU:Z:${OLD_RGs[1]}\b//g; "
          fi
          break
        fi
      done
      if [[ $FOUND -eq 0 ]]; then
        echo "Error: actual ID: $ACTUAL_ID for BAM RG ID: ${OLD_RGs[0]} cannot find a match in the manifest: $MANIFEST!"
        exit 1
      fi
   fi
   
done
echo
SED_CMD=$SED_CMD1" "$SED_CMD2" s/\tLB:Z:N\/A//g;"
CMD="samtools view -h $IN_BAM | sed '$SED_CMD' | samtools view -S -b - > ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: samtools is failed to modify the body part!"
   exit 1
fi



CMD="samtools view -H ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam | awk -v ids=\"$RG\" -v sm=\"$SM\" 'BEGIN{first=1}{if (\$1~/^@RG/) {
                                                if (first==1){
                                                split(ids,aa,\";\"); 
                                                for (i in aa){
                                                  split(aa[i],bb,\"|\");
                                                  print \"@RG\tID:\"bb[1]\"\tSM:\"sm\"\tLB:\"bb[3]\"\tPL:ILLUMINA\tPU:\"bb[2]\"\tCN:CGR\"
                                                }
                                                first=0
                                                }} 
                                                else{
                                                  print 
                                                 }}' > ${TMP_DIR}/${ANALYSIS_ID}_header.sam"
echo
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: samtools is failed to modify the header part!"
   exit 1
fi



CMD="samtools reheader ${TMP_DIR}/${ANALYSIS_ID}_header.sam ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam > ${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader.bam"
echo
echo $CMD
eval $CMD

if [[ $? -ne 0 ]]; then
   echo "Error: samtools reheader is failed!"
   exit 1
fi

echo
CMD="rm -f ${TMP_DIR}/${ANALYSIS_ID}_tmp.bam ${TMP_DIR}/${ANALYSIS_ID}_header.sam "
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: removing the temp files is failed!"
   exit 1
fi

echo
CMD="samtools sort ${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader.bam -o ${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader_sorted.bam"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: samtools sort is failed!"
   exit 1
fi

echo
CMD="samtools index ${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader_sorted.bam "
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: samtools index is failed!"
   exit 1
fi

CMD="rm -f ${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader.bam"
echo $CMD
eval $CMD

echo
echo "Comparing the counts between the In and Out BAMs.."
OLD_COUNT=`samtools view -c $IN_BAM`
NEW_COUNT=`samtools view -c ${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader_sorted.bam`
echo $OLD_COUNT" "$NEW_COUNT
if [[ $OLD_COUNT -ne $NEW_COUNT ]]; then
   echo "Error: the count is inconsistent!"
   exit 1
fi

echo "[$(date)]: Remove duplicates for '$ANALYSIS_ID'..."
CMD="java -Xmx8g -jar $PICARD MarkDuplicates \
        INPUT=${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader_sorted.bam \
        OUTPUT=$OUT_BAM \
        METRICS_FILE=$OUT_DEDUP_metrics
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT \
        CREATE_INDEX=true"
echo $CMD
eval $CMD

if [[ $? -ne 0 ]]; then
    echo "Error: $(date) dedup is failed!"
    exit 1
fi

CMD="rm -f ${TMP_DIR}/${ANALYSIS_ID}_tmp_reheader_sorted.bam"
echo $CMD
eval $CMD


END=`date +%s`

RUNTIME=$((END-START))
echo "Runtime: $RUNTIME"
echo "[$(date)]: done!"
touch ${LOG_DIR}/fix_${ANALYSIS_ID}.DONE
rm -f ${LOG_DIR}/fix_${ANALYSIS_ID}.WORKING
