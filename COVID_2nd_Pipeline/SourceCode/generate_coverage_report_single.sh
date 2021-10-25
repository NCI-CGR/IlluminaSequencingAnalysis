#!/bin/sh
set -o pipefail
# source global configuration files
SCRIPT=$(readlink -f "$0")
DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc

# takes in arguments
IN_BAM=$1
OUT_REPORT=$2
MANIFEST=$3
strFlagWorking=$4
strFlagDone=$5
#MERGE_LOG=$4

echo ====================
date
echo Coverage Report Generation
echo IN_BAM: $IN_BAM
echo Append Stats To: $OUT_REPORT
echo ====================
if [[ $IN_BAM == *bqsr_final_out* ]]; then
IN_BAM_BASE=`basename $IN_BAM _bqsr_final_out.bam | cut -d_ -f1-`
ANALYSIS_ID=`basename $IN_BAM _bqsr_final_out.bam | cut -d_ -f2-`
else
IN_BAM_BASE=`basename $IN_BAM .bam | cut -d_ -f1-`
ANALYSIS_ID=`basename $IN_BAM .bam | cut -d_ -f2-`
fi
#kevin's glu script will pend a NOID to the sample name if the fields is not same as familila sample naming convention. Get rid of the NOID to match back the manifest file.
ANALYSIS_ID=`echo "${ANALYSIS_ID//_NOID}"`
#retrieve limsID from manifest file
#LIMSID=`awk -F '\t' -v QUERY_SAMPLE=$ANALYSIS_ID '{if (QUERY_SAMPLE ~ $9) {print $8}}' $MANIFEST`
LIMSID=`awk -F ',' -v QUERY_SAMPLE=$ANALYSIS_ID '{if (QUERY_SAMPLE ~ $13) {print $8}}' $MANIFEST`
LIMSID=`echo $LIMSID | cut -d" " -f1`


#SR_SUBJECT_ID=`awk -F '\t' -v QUERY_SAMPLE=$ANALYSIS_ID '{if ($9 == QUERY_SAMPLE) {print $11}}' $MANIFEST | sort | uniq`
SR_SUBJECT_ID=`awk -F ',' -v QUERY_SAMPLE=$ANALYSIS_ID '{if ($13 == QUERY_SAMPLE) {print $14}}' $MANIFEST | sort | uniq`
SR_SUBJECT_ID=`echo $SR_SUBJECT_ID | sed 's/ /,/g'|sed 's/\r//g'`
echo $LIMSID,$SR_SUBJECT_ID
REPORT_LINE="$(date +"%m/%d/%Y")\t${LIMSID}\t${SR_SUBJECT_ID}\t${IN_BAM_BASE}"


#retrieve captureKit from manifest file and use the corresponding bed file
#CAPTUREKIT=`awk -F '\t' -v QUERY_SAMPLE=$LIMSID -v ANALYSISID=$ANALYSIS_ID '{if ($8 == QUERY_SAMPLE && $9 == ANALYSISID )  {print $10}}' $MANIFEST|sort|uniq`
CAPTUREKIT=`awk -F ',' -v ANALYSISID=$ANALYSIS_ID '{if ($13 == ANALYSISID )  {print $12}}' $MANIFEST|sort|uniq`
CAPTUREKIT=`echo $CAPTUREKIT|sed 's/ //g'`


echo "CAPTUREKIT: ${CAPTUREKIT}"
if [[ $CAPTUREKIT == *3.0* ]]; then
		EXOME_TARGETS=$EXOME_TARGETS_v3
		EXOME_TARGETS_TOTAL_BASES=$EXOME_TARGETS_TOTAL_BASES_v3
elif [[ $CAPTUREKIT == *UTR* ]]; then
		EXOME_TARGETS=$EXOME_TARGETS_v3plusUTR
		EXOME_TARGETS_TOTAL_BASES=$EXOME_TARGETS_TOTAL_BASES_v3plusUTR
elif [[ $CAPTUREKIT == *EZ_Choice_TL-Kid-Lung* ]]; then
		EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/161006_HG19_Kid-Lung-Extra_EZ_HX3_primary_targets.bed	
		EXOME_TARGETS_TOTAL_BASES=1780492
elif [[ $CAPTUREKIT == EZ_Choice_Taylor-ESCC-Kid-Lung-Extra ]]; then
		EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/EZ_Choice_Taylor-ESCC-Kid-Lung-Extra.bed
		EXOME_TARGETS_TOTAL_BASES=1939177
elif [[ $CAPTUREKIT == EZ_Choice_Kid-Lung-Extra ]]; then
		EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/Kid-Lung-Extra_capture_targets.bed
		EXOME_TARGETS_TOTAL_BASES=1897580
elif [[ $CAPTUREKIT == EZ_Choice_CGB-368genes ]]; then
		EXOME_TARGETS=/DCEG/CGF/Laboratory/LIMS/Oligo_Orders/NimblegenCapture/CGB-368/CGB-368genes_primary_targets.bed
		EXOME_TARGETS_TOTAL_BASES=2674426
elif [[ $CAPTUREKIT == EZ_Choice_Dean-Koshiol-4 ]]; then
		EXOME_TARGETS=/DCEG/CGF/Laboratory/LIMS/Oligo_Orders/NimblegenCapture/Dean-Koshiol-4/Selection_Results/Dean-Koshiol-4_capture_targets.bed
		EXOME_TARGETS_TOTAL_BASES=2956332
elif [[ $CAPTUREKIT == EZ_Choice_Karina-XP ]]; then
    EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/Karina-XP_primary_targets.bed
		EXOME_TARGETS_TOTAL_BASES=941427
elif [[ $CAPTUREKIT == EZ_Choice_Dean-MexPC ]]; then
    EXOME_TARGETS=/DCEG/CGF/Laboratory/LIMS/Oligo_Orders/NimblegenCapture/Dean-MexPC/Selection_Results/Dean-MexPC_primary_targets.bed
		EXOME_TARGETS_TOTAL_BASES=2252101
elif [[ $CAPTUREKIT == *_WGS_* ]]; then
    EXOME_TARGETS=${WGS_BED}
		EXOME_TARGETS_TOTAL_BASES=${WGS_TOTAL_BASES}
else
	echo no bed file found!
	exit 1
#	EXOME_TARGETS=$EXOME_TARGETS_v3
#	EXOME_TARGETS_TOTAL_BASES=$EXOME_TARGETS_TOTAL_BASES_v3
fi

# EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/variant_scripts_V3_in_process/variant_intervals_beds/variant_calling_intervals_NV3UTR_NV3_250padded_corrected_4000parts/NV3UTR_NV3_250padded_corrected.bed
# EXOME_TARGETS_TOTAL_BASES=208400278

echo "EXOME_TARGETS            : ${EXOME_TARGETS}"
echo "EXOME_TARGETS_TOTAL_BASES: ${EXOME_TARGETS_TOTAL_BASES}"

for REFERENCE in $REFERENCE_CDS $EXOME_TARGETS ; do
#for REFERENCE in /DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/150911_HG19_DS-Neuro_EZ_HX3_capture_targets.bed;do
	if [[ $REFERENCE == $EXOME_TARGETS ]]; then
		TOTAL_BASES=$EXOME_TARGETS_TOTAL_BASES
	else
		TOTAL_BASES=$REFERENCE_CDS_TOTAL_BASES
	fi

	echo "Reference BED File: $REFERENCE; Total Bases: $TOTAL_BASES"

	# generate the base-by-base pileup based on the designated bed file, extract the most essential
	# information (chromosome, position, read depth) into the report file
	echo "Building pileup for coverage analysis..."

	TMP_FILE=$SCRATCH_DIR/${IN_BAM_BASE}.coverage
	echo "Temp File: $TMP_FILE"
	rm -rf $TMP_FILE
	####
        #use option -A if want to included filtered orphan reads
	samtools depth -@ 8 -d 0 -q 0 -b $REFERENCE $IN_BAM  > $TMP_FILE
	####

	if [ ! -f $TMP_FILE ]; then
		echo "Output report from samtools mpileup $TMP_FILE expected but not found." 1>&2
		exit 1
	fi

	for THRESHOLD in 1 5 10 15 50; do
		CMD="awk -v threshold=$THRESHOLD '{if(\$3>=threshold) L++} END {print L}' $TMP_FILE"
		#echo $CMD
		TOTAL_BASES_GT=`eval $CMD`
		PERCENT_BASES_GT=`echo $TOTAL_BASES_GT/$TOTAL_BASES | bc -l  | xargs printf "%5.4f" `
		#AVG_COVERAGE=`awk '{sum+=$3} END {print sum/NR}' $REPORT_FILE | xargs printf "%2.2f" `
		echo "Total Exome Bases >= ${THRESHOLD}x Coverage: $TOTAL_BASES_GT; Percent: $PERCENT_BASES_GT"

		REPORT_LINE="${REPORT_LINE}\t${TOTAL_BASES_GT}\t${PERCENT_BASES_GT}"
	done
	AVG_COVERAGE=`awk -v total=$TOTAL_BASES '{sum+=$3} END {print sum/total}' $TMP_FILE | xargs printf "%2.2f" `
  REPORT_LINE="${REPORT_LINE}\t${AVG_COVERAGE}"
done
# TOTAL_READS=`grep "pairs never matched" $MERGE_LOG |cut -d" " -f3`

# DUPLICATE_READS=`grep "records as duplicates" $MERGE_LOG |cut -d" " -f3`
# PERCENT_DUP=`echo $DUPLICATE_READS/$TOTAL_READS |bc -l |xargs printf "%5.4f"`
# REPORT_LINE="${REPORT_LINE}\t${PERCENT_DUP}"

# OPTICAL_DUPLICATE_READS=`grep "optical duplicate clusters" $MERGE_LOG |cut -d" " -f3`
# PERCENT_OPTICAL_DUPLICATE_READS=`echo $OPTICAL_DUPLICATE_READS/$TOTAL_READS |bc -l |xargs printf "%5.4f"`
# REPORT_LINE="${REPORT_LINE}\t${PERCENT_OPTICAL_DUPLICATE_READS}"

# AVG_COVERAGE=`awk -v total=$TOTAL_BASES '{sum+=$3} END {print sum/total}' $TMP_FILE | xargs printf "%2.2f" `
# REPORT_LINE="${REPORT_LINE}\t${AVG_COVERAGE}"
echo -e $REPORT_LINE
rm -rf $TMP_FILE

# now print out all stats to the report file
#echo -e $REPORT_LINE  >> $OUT_REPORT

echo ====================
date
echo Done!
echo ====================

# update flag
rm ${strFlagWorking}
touch ${strFlagDone}
