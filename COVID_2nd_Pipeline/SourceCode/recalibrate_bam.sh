#!/bin/sh
# . ${SCRIPT_DIR:-.}/global_config_bash.rc
. ${SCRIPT_DIR}/global_config_bash.rc
module load jdk/1.8.0_111
module load samtools/1.8
IN_BAM=$1
OUT_BAM=$2
VARIANT_TYPE=$3

IN_DIR=`dirname $IN_BAM`
IN_NAME=`basename $IN_BAM .bam`
# the input index could be name.bai or name.bam.bai
# try to guess it slightly smartly
if [[ -f ${IN_DIR}/${IN_NAME}.bai ]] ; then
    IN_BAI=${IN_DIR}/${IN_NAME}.bai
else
    IN_BAI=$IN_BAM.bai
fi

OUT_DIR=`dirname $OUT_BAM`
OUT_NAME=`basename $OUT_BAM .bam`
OUT_BAI=${OUT_DIR}/${OUT_NAME}.bam.bai
if [[ ! -d $OUT_DIR ]]; then
  mkdir -p $OUT_DIR
fi

if [[ $? -ne 0 ]]; then
   echo "Error: mkdir $OUT_DIR is failed!"
   exit 1
fi
# FLG_WORKING=$OUT_BAM.$RECALIBRATION_WORKING_FLAG_EXTENSION
# FLG_INQUEUE=$OUT_BAM.$RECALIBRATION_INQUEUE_FLAG_EXTENSION

FLG_WORKING=${BAM_REFORMATTED_RECALIBRATED_DIR}/${IN_NAME}.bam.${RECALIBRATION_WORKING_FLAG_EXTENSION}
FLG_INQUEUE=${BAM_REFORMATTED_RECALIBRATED_DIR}/${IN_NAME}.bam.${RECALIBRATION_INQUEUE_FLAG_EXTENSION}

touch $FLG_WORKING
rm -f $FLG_INQUEUE


LEFT_ALIGNED_BAM=${TMP_DIR}/${OUT_NAME}_left_aligned.bam
LEFT_ALIGNED_BAI=${TMP_DIR}/${OUT_NAME}_left_aligned.bai

REALIGN_TARGETS=${TMP_DIR}/${OUT_NAME}_realign.list
REALIGN_BAM=${TMP_DIR}/${OUT_NAME}_realign.bam
REALIGN_BAI=${TMP_DIR}/${OUT_NAME}_realign.bai

BQSR_TARGETS=${TMP_DIR}/${OUT_NAME}_recal.list
BQSR_BAM=${TMP_DIR}/${OUT_NAME}_bqsr.bam
BQSR_BAI=${TMP_DIR}/${OUT_NAME}_bqsr.bai

OUT_TMP_BAM_PREFIX=${TMP_DIR}/${OUT_NAME}_final_out
OUT_TMP_BAM=$OUT_TMP_BAM_PREFIX.bam
OUT_TMP_BAI=$OUT_TMP_BAM_PREFIX.bam.bai

OUT_BAM_PREFIX_SOMATIC=${TMP_DIR}/${OUT_NAME}_somatic_final_out
OUT_BAM_SOMATIC=${OUT_DIR}/${OUT_NAME}_bqsr_final_out.bam
OUT_BAI_SOMATIC=${OUT_DIR}/${OUT_NAME}_bqsr_final_out.bam.bai


echo "Clean up the potential leftover ...."
CMD="rm -f ${TMP_DIR}/${OUT_NAME}* ${TMP_DIR}/${OUT_NAME}*.* ${TMP_DIR}/${OUT_NAME}*.*.*"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
    echo "Error: rm is failed!"
    exit 1
fi


echo ========================
date
echo Recalibrate BAM according to reference InDels and SNPs
echo Input: $IN_BAM
echo Output: $OUT_BAM
echo Type: $VARIANT_TYPE
echo ========================

echo
# Make sure all InDels are left-aligned
# turned on by default
if [[ $LEFT_ALIGN_INDELS != "false" ]] ; then
        echo "[$(date)] Left align InDels..."
        CMD="java -Xmx8g -jar $GATK   \
                -T LeftAlignIndels   \
                -R $REFERENCE_GENOME \
                -I $IN_BAM \
                -o $LEFT_ALIGNED_BAM"
	echo $CMD
	eval $CMD
	if [[ $? -ne 0 ]]; then
	    echo "Error: LeftAlignIdels is failed!"
	    exit 1
	fi
        if [[ ! -f $LEFT_ALIGNED_BAM ]] ; then
            echo "Error: Expected output $LEFT_ALIGNED_BAM is missing."
            exit 1
        fi

        echo "[$(date)] Left align is done."

else
        echo "[$(date)] InDel left alignment skipped as LEFT_ALIGN_INDELS is set to false. Input BAM file assumed to be left-aligned."
        rm -f $LEFT_ALIGNED_BAM $LEFT_ALIGNED_BAI
        ln -s $IN_BAM $LEFT_ALIGNED_BAM
fi

echo
echo "[$(date)] Calculating target intervals that need local realignment..."
CMD="java -Xmx8g -jar $GATK        \
	-T RealignerTargetCreator \
	-R $REFERENCE_GENOME      \
	-known $REFERENCE_INDELS  \
	-known $REFERENCE_MILLS_INDELS  \
	-I $LEFT_ALIGNED_BAM      \
	-o $REALIGN_TARGETS \
        --disable_auto_index_creation_and_locking_when_reading_rods"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
     echo "Error: GATK RealignerTargetCreator is failed!"
     exit 1
fi
if [[ ! -f $REALIGN_TARGETS ]] ; then
     echo "Error: Expected output $REALIGN_TARGETS is missing." 1>&2
     exit 1
fi
echo "[$(date)] Calculating target intervals is done!"
echo
echo "[$(date)] Carrying out local alignment in target intervals..." &&
CMD="java -Xmx8g -jar $GATK                 \
	-T IndelRealigner                  \
	-R $REFERENCE_GENOME               \
	-targetIntervals $REALIGN_TARGETS  \
	-known $REFERENCE_INDELS           \
	-known $REFERENCE_MILLS_INDELS     \
        -rf NotPrimaryAlignment            \
	-I $LEFT_ALIGNED_BAM               \
	-o $REALIGN_BAM"
	
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: GATK IndelRealigner is failed!"
   exit 1
fi
if [[ ! -f $REALIGN_BAM ]] ; then
   echo "Error: Expected output $REALIGN_BAM is missing." 1>&2
   exit 1
fi

#HAS_DONE=`tail -n 10 ${CLUSTER_JOB_LOG_DIR}/_recalibration_${IN_NAME}.bam.stderr | grep "Total runtime" | wc -l`
HAS_DONE=`tail -n 10 ${CLUSTER_JOB_LOG_DIR}/*/_recalibration_${IN_NAME}.bam.stderr | grep "Total runtime" | wc -l`
if [[ $HAS_DONE -eq 0 ]] ; then
     echo "Error: Expected output $REALIGN_BAM is truncated." 1>&2
     exit 1
fi

echo "[$(date)] Local realignment is done!"

# delete old version of staging file to save space
CMD="rm -f $IN_BAM_WITH_READGROUP $IN_BAI_WITH_READGROUP $REORDERED_BAM $REORDERED_BAI $LEFT_ALIGNED_BAM $LEFT_ALIGNED_BAI $REALIGN_TARGETS"
echo $CMD
eval $CMD


echo
echo "[$(date)] Doing the properly paired filtering ..."
CMD="samtools view -b -f 3 $REALIGN_BAM > $OUT_TMP_BAM"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
    echo "Error: samtools filtering is failed!"
    exit 1
fi
echo "[$(date)] properly-paried filtering is done!"
echo
CMD="samtools sort -T ${OUT_TMP_BAM_PREFIX}_sorted_tmp -o ${OUT_TMP_BAM_PREFIX}_sorted.bam $OUT_TMP_BAM"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
     echo "Error: samtools sorting is failed!"
     exit 1
fi
echo
echo "[$(date)] Doing indexing ..."
CMD="samtools index ${OUT_TMP_BAM_PREFIX}_sorted.bam ${OUT_TMP_BAM_PREFIX}_sorted.bai" 
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
     echo "Error: samtools indexing is failed!"
     exit 1
fi
echo "[$(date)] Indexing is done!"
echo

CMD="mv ${OUT_TMP_BAM_PREFIX}_sorted.bam $OUT_BAM"
echo $CMD
eval $CMD

if [[ $? -ne 0 ]]; then
     echo "Error: mv the BAM file is failed!"
     exit 1
fi

CMD="mv ${OUT_TMP_BAM_PREFIX}_sorted.bai $OUT_BAI"
echo $CMD
eval $CMD

if [[ $? -ne 0 ]]; then
     echo "Error: mv the BAI file is failed!"
     exit 1
fi
      
CMD="rm -f ${OUT_TMP_BAM_PREFIX}_sorted_tmp* $OUT_TMP_BAM $OUT_TMP_BAI"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
    echo "Error: removing the temp files is failed!"
    exit 1
fi



# by BQSR by default turned off, as it takes too long time
if [[ $VARIANT_TYPE == "SOMATIC" ]]; then
        echo
	echo "[$(date)] This BAM is for SOMATIC project. Calculating target sites that need base quality score recalibration (BQSR)..."
	CMD="java -Xmx16g -jar $GATK               \
		-T BaseRecalibrator              \
		-R $REFERENCE_GENOME             \
		-knownSites $REFERENCE_SNPS2      \
		-knownSites $REFERENCE_MILLS_INDELS      \
		-knownSites $REFERENCE_INDELS    \
		-I $REALIGN_BAM                  \
		-o $BQSR_TARGETS"
        echo $CMD
        eval $CMD

	if [[ $? -ne 0 ]]; then
	  echo "Error: failed in BQSR!"
	  exit 1
	fi
        echo
	echo "[$(date)] BQSR Target sites calculation is done."

	echo "[$(date)] Carrying out base quality score recalibration (BQSR)..."
	CMD="java  -Xmx16g -jar $GATK             \
		-T PrintReads                   \
		-R $REFERENCE_GENOME            \
		-BQSR $BQSR_TARGETS             \
		-I $REALIGN_BAM                 \
		-o $BQSR_BAM"
        echo $CMD
        eval $CMD

	if [[ $? -ne 0 ]]; then
	  echo "Error: failed in BQSR!"
	  exit 1
	fi

	echo "[$(date)] BQSR is Done."
	
	# BQSR unfortunately generate query-name-sorted BAM
	# Convert the realigned query-name-sorted BAM to coordinate-sorted bams
        echo
	echo "[$(date)] Converting BAM to coordinate-sorted format..."
           
        
        CMD="samtools sort -T ${OUT_BAM_PREFIX_SOMATIC}_sorted_tmp -o ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam $BQSR_BAM" 
        echo $CMD
        eval $CMD
        if [[ $? -ne 0 ]]; then
            echo "Error: samtools sort is failed!"
            exit 1
        fi
	if [[ ! -f ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam ]] ; then
	    echo "Error: Expected output ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam is missing."
	    exit 1
	fi
	
        CMD="samtools index ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam.bai"
        echo $CMD
        eval $CMD
        if [[ $? -ne 0 ]]; then
            echo "Error: samtools index is failed!"
            exit 1
        fi
        if [[ ! -f ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam.bai ]] ; then
            echo "Error: Expected output ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam.bai is missing."
            exit 1
        fi

        echo 
	CMD="rm -rf $BQSR_TARGETS $BQSR_BAM $BQSR_BAI ${OUT_BAM_PREFIX_SOMATIC}_sorted_tmp*"
        echo $CMD
        eval $CMD

        
        echo

        CMD="mv ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam $OUT_BAM_SOMATIC"
        echo $CMD
        eval $CMD
        echo
        CMD="mv ${OUT_BAM_PREFIX_SOMATIC}_sorted.bam.bai $OUT_BAI_SOMATIC" 
        echo $CMD
        eval $CMD

        echo "[$(date)] doing the last piece of housekeeping for $OUT_BAM_SOMATIC ..."
        CMD="rm -f ${OUT_BAM_PREFIX_SOMATIC}_sorted*"
        echo $CMD
        eval $CMD     
        chgrp ncicgf_dceg_exome $OUT_BAM_SOMATIC 2>/dev/null && chmod g+r $OUT_BAM_SOMATIC
        chgrp ncicgf_dceg_exome $OUT_BAI_SOMATIC 2>/dev/null && chmod g+r $OUT_BAI_SOMATIC
        echo "[$(date)] $OUT_BAM_SOMATIC is done!"

fi


echo
echo "[$(date)] Doing the last piece of housekeeping..." 
chgrp ncicgf_dceg_exome $OUT_BAM && chmod g+r $OUT_BAM 2>/dev/null
chgrp ncicgf_dceg_exome $OUT_BAI && chmod g+r $OUT_BAI 2>/dev/null
CMD="rm -f $REALIGN_BAM $REALIGN_BAI"
echo $CMD
eval $CMD
CMD="rm -f $FLG_WORKING"
echo $CMD
eval $CMD

echo "[$(date)] All done!"

