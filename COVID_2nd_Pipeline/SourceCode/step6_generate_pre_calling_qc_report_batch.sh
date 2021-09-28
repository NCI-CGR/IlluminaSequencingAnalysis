#!/bin/sh

#SCRIPT=$(readlink -f "$0")
#DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
#. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc
. ./global_config_bash.rc

DATE=`echo $(date +%m%d%Y)`
LOG_DIR=${CLUSTER_JOB_LOG_DIR}/${DATE}
mkdir -p $LOG_DIR 2>/dev/null
rm ./nohup100.out

MANIFEST_FILE=$1
shift

if [[ $MANIFEST_FILE == "" ]];then
	echo "Please provide manifest file!"
	exit 0
fi

#test if the manifest file is txt file, change to txt if the file is csv
if [[ ${MANIFEST_FILE: -4} == ".csv" ]];then
	sed -e 's/,/\t/g' $MANIFEST_FILE >${BUFFER_DIR:-.}/pre_calling_qc_manifest_temp.txt
	MANIFEST_FILE=${BUFFER_DIR:-.}/pre_calling_qc_manifest_temp.txt
fi

# Some of the batch reporting jobs may not yield result
# as expected (e.g. due to dead cluster node).
# In such case, we'll want to rerun this script, and use
# the last "good" report file as an argument following the script 
# (i.e. "sh ./generate_coverage_report_batch.sh <my_last_report.txt>"
# The script will then enter a "patch mode", which essentially search for
# the IDs that are NOT in <my_last_report.txt> and then submit the cluster job accordingly.
# All IDs that are already in <my_last_report.txt> will then be skipped

if [[ $# == 1 ]]; then
	# if there is an argument, enter patch mode
	PATCH_MODE=true
	PREQC_REPORT_FILE=$1
	echo "Entering patch mode to fill up the \"holes\" in existing coverage report file $PREQC_REPORT_FILE"
else
	# otherwise automatically generate the report file name
	PATCH_MODE=false
	PREQC_REPORT_FILE=${COVERAGE_REPORT_DIR}/pre_calling_qc_report_${DATE}.txt
	
	REPORT_LINE="Analysis ID\tAVERAGE_READ_LENGTH\tPERCENT_TOTAL_READS_WITH_3LESS_CUT\tPF_HQ_MEDIAN_MISMATCHES\tPF_HQ_MISMATCH_RATE_F\tPF_HQ_MISMATCH_RATE_R\tSTRAND_BALANCE_F\tSTRAND_BALANCE_R\tBAD_CYCLES\tPCT_CHIMERAS\tPCT_ADAPTER"\
"\tPERCENT_PF_HQ_ALIGNED_READS\tPERCENT_PF_HQ_ALIGNED_BASES_L\tPERCENT_PF_HQ_ALIGNED_Q20_BASES_L\tPERCENT_PF_HQ_ALIGNED_BASES_R\tPERCENT_PF_HQ_ALIGNED_Q20_BASES_R\tBASES_QUAL_AVE"\
"\tAT_DROPOUT\tGC_DROPOUT\tHS_AT_DROPOUT\tHS_GC_DROPOUT\tFOLD_80_BASE_PENALTY\tMEDIAN_INSERT_SIZES\tMEAN_INSERT_SIZES\tAVG_COVERAGE"\
"\tLOWEST_PREADATPER_TOTAL_QSCORE_BASE\tLOWEST_PREADATPER_TOTAL_QSCORE\tLOWEST_BAITBIAS_TOTAL_QSCORE_BASE\tLOWEST_BAITBIAS_TOTAL_QSCORE\tOXIDATION_Q_AVE\tLOWEST_OXIDATION_Q\tLOWEST_OXIDATION_Q_CONTEXT"\
"\tPERCENT_LOW_COVERAGE_CAPTURE_LOCI\tPERCENT_NO_COVERAGE_CAPTURE_LOCI\tPERCENT_POOR_MAPPING_QUALITY_LOCI\tPERCENT_READS_ENDING_IN_INDELS\tPERCENT_READS_ENDING_IN_SOFTCLIPS"

    echo -e $REPORT_LINE > $PREQC_REPORT_FILE
fi	
	

# all new BAM files are generated under the incoming directory
# these are typically the files that need coverage report
echo "awk -F "\t" -v col='GROUP' '{for (i=1;i<=NF;i++) if ($i==col) print i; exit}' $MANIFEST_FILE"
GROUP_COL=$(awk -F "\t" -v col='GROUP' '{for (i=1;i<=NF;i++) if ($i==col) print i; exit}' $MANIFEST_FILE)
echo $GROUP_COL
ANALYSISID_COL=$(awk -F "\t" -v col='ANALYSIS ID' '{for (i=1;i<=NF;i++) if ($i~col) print i; exit}' $MANIFEST_FILE)
echo $ANALYSISID_COL
BAM_DIR=$(dirname ${MANIFEST_FILE})/../bam_location
cd $BAM_DIR
for BAM in $(awk -F "\t" -v group=$GROUP_COL -v analysisid=$ANALYSISID_COL '{if(NR>1){print $group"/"$group"_"$analysisid".bam"}}' $MANIFEST_FILE | sort | uniq); do	
		echo $BAM
#for BAM in $(cat GIAB.lst); do
# in case we may have forgotten to generate this pre-calling qc report
# before we move incoming BAM files to original BAM folder
# here is another way, we'll generate the reports for all BAM files

#for BAM in ${BAM_DIR}/*/*.bam; do

# for benchmarking purposeSH
#for BAM in /DCEG/Projects/Exome/SequencingData/BAM_reduced_coverage_test/*.bam; do
	#echo "==> $BAM <=="
	NAME=`basename $BAM .bam`
#	IN_BAM=/DCEG/Projects/Exome/SequencingData/BAM_recalibrated_dedup/*/${NAME}.bam
	#echo $NAME
#	BAM=${BAM_INCOMING_DIR}/../BAM_original/OSTE/${NAME}.bam
	if [[ $PATCH_MODE == "true" ]] ; then
		if [[ $(awk -F "\t" -v searchid="$NAME" '{if ($1 == searchid) {print} }' $PREQC_REPORT_FILE | wc -l) == 0 ]]; then
			DO_THIS_BAM=true
			echo "$NAME not in existing report file. Redoing the report generation..."
		else
			DO_THIS_BAM=false
			#echo "$NAME already in report file. Skipped."
		fi
	else 
		DO_THIS_BAM=true
	fi

	if [[ ($DO_THIS_BAM == "true") || ($PATCH_MODE == "false") ]] ; then
		CMD="qsub -q $QUEUE -cwd \
			-o ${LOG_DIR}/_pre_calling_qc_report_${NAME}.stdout -e ${LOG_DIR}/_pre_calling_qc_report_${NAME}.stderr \
			-N PRECALLING_QC.$NAME \
			-S /bin/sh ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/pre_calling_qc_single.sh $BAM $PREQC_REPORT_FILE $MANIFEST_FILE"
		echo $CMD
		
		echo $CMD >> ./nohup100.out
		eval $CMD
	fi
done
