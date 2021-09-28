#!/bin/sh


if [[ $# -lt 2 ]]; then
  echo "Error: please run it with a least 2 argument (Manifest file, and build_BAM dir)!"
  exit 1
fi

SCRIPT=$(readlink -f "$0")
DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc

DATE=`echo $(date +%m%d%Y)`
LOG_DIR=${CLUSTER_JOB_LOG_DIR}/${DATE}

# Get two most important arguments -->
#MANIFEST=/DCEG/Projects/Exome/builds/build_SR0407-004_Data_Delivery_2019_23058/Manifest/NP0407-HE9-ANALYSIS-MANIFEST.csv
MANIFEST=$1
dirBuildBAM=$2
#create BAM list file
strBAMList=${dirBuildBAM}/bam.lst
if [ -f "${strBAMList}" ]; then
  rm ${strBAMList}
fi
# Create bam.list
find ${dirBuildBAM} -iname '*.bam' > ${strBAMList}
#<--
#MERGE_LOG_DIR=/DCEG/Projects/Exome/SequencingData/variant_scripts/logs/GATK/patch_build_bam_20190321
mkdir -p $LOG_DIR 2>/dev/null

# Some of the batch reporting jobs may not yield result
# as expected (e.g. due to dead cluster node).
# In such case, we'll want to rerun this script, and use
# the last "good" report file as an argument following the script
# (i.e. "sh ./generate_coverage_report_batch.sh <my_last_report.txt>"
# The script will then enter a "patch mode", which essentially search for
# the IDs that are NOT in <my_last_report.txt> and then submit the cluster job accordingly.
# All IDs that are already in <my_last_report.txt> will then be skipped

if [[ $# == 3 ]]; then
	# if there is an argument, enter patch mode
	PATCH_MODE=true
	COVERAGE_REPORT_FILE=$3
	echo "Entering patch mode to fill up the \"holes\" in existing coverage report file $COVERAGE_REPORT_FILE"
else
	# otherwise automatically generate the report file name
	PATCH_MODE=false

	#Init folder -->
	if [ ! -d ${COVERAGE_REPORT_DIR} ]; then
	  mkdir -p ${COVERAGE_REPORT_DIR}
	fi
  #<--

	COVERAGE_REPORT_FILE=${COVERAGE_REPORT_DIR}/coverage_report_${DATE}.txt

	# output report header, note this needs to match the actual numbers generated in ./generate_coverage_report_single.sh
	 REPORT_LINE="Report Date\tLIMS ID\tSR SUBJECT ID\tAnalysis ID" > $COVERAGE_REPORT_FILE
	 for REFERENCE in UCSC CaptureKit; do
		 for THRESHOLD in 1 5 10 15 50; do
	#REPORT_LINE="Report Date\tAnalysis ID" > $COVERAGE_REPORT_FILE
	#for REFERENCE in secondary_analysis_bed_files/*.bed; do
	#	for THRESHOLD in 15; do
	#		REFERENCE2=$(basename $REFERENCE .bed)
			REPORT_LINE="${REPORT_LINE}\tBases ${REFERENCE} >= ${THRESHOLD}X\t% ${REFERENCE} >= ${THRESHOLD}X"
	#		REPORT_LINE="${REPORT_LINE}\t% ${REFERENCE2} >= ${THRESHOLD}X\t ${REFERENCE2} Avarage Coverage"
		  done
	  done

	#REPORT_LINE="${REPORT_LINE}\tUCSC Avarage Coverage"
	REPORT_LINE="${REPORT_LINE}\t% Merge Dup"
	echo -e $REPORT_LINE > $COVERAGE_REPORT_FILE
	REPORT_LINE="${REPORT_LINE}\t% Merge Optical Dup"
	echo -e $REPORT_LINE > $COVERAGE_REPORT_FILE
	REPORT_LINE="${REPORT_LINE}\tCaptureKit Avarage Coverage"
	echo -e $REPORT_LINE > $COVERAGE_REPORT_FILE
fi

# all new BAM files are generated under the incoming directory
# these are typically the files that need coverage report
#for BAM in $BAM_INCOMING_DIR/*/*.bam; do

#for BAM in `cat /DCEG/Projects/Exome/SequencingData/sync_script_v3/LB_deduplication_test/recalibrated/ec.lst`; do
#for BAM in `cat /DCEG/Projects/Exome/builds/build_SR0407-004_Data_Delivery_2019_23058/bam.lst`; do
for BAM in `cat ${strBAMList}`; do
	NAME=`basename $BAM .bam`
	if [[ $BAM == *bqsr_final_out* ]]; then
  	BASE_NAME=`basename $BAM _bqsr_final_out.bam `
#  	MERGE_LOG=${MERGE_LOG_DIR}/_build_bam_${BASE_NAME}.stderr
#	else
#  	MERGE_LOG=${MERGE_LOG_DIR}/_build_bam_${NAME}.stderr
	fi
	
#	IN_BAM=/DCEG/Projects/Exome/SequencingData/BAM_recalibrated_dedup/*/${NAME}.bam
	#echo $NAME
#	BAM=${BAM_INCOMING_DIR}/../BAM_original/OSTE/${NAME}.bam
	if [[ $PATCH_MODE == "true" ]] ; then
		if [[ $(awk -F "\t" -v searchid="$NAME" '{if ($2 == searchid) {print} }' $COVERAGE_REPORT_FILE | wc -l) == 0 ]]; then
			DO_THIS_BAM=true
			echo "$NAME not in existing report file. Redoing the report generation..."
		else
			DO_THIS_BAM=false
			#echo "$NAME already in report file. Skipped."
		fi
	fi

  if [[ ($DO_THIS_BAM == "true") || ($PATCH_MODE == "false") ]] ; then
#		CMD="cd $DCEG_SEQ_POOL_SCRIPT_DIR; qsub -q $QUEUE -cwd \
#			      -o ${LOG_DIR}/_coverage_report_${NAME}.stdout -e ${LOG_DIR}/_coverage_report_${NAME}.stderr \
#			      -N COVREP.$NAME \
#			      -S /bin/sh ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/generate_coverage_report_single.sh $BAM $COVERAGE_REPORT_FILE $MANIFEST" #$MERGE_LOG "
    CMD="cd $DCEG_SEQ_POOL_SCRIPT_DIR; \
         sbatch --ntasks=1 \
                --nodes=1 \
                --job-name=COVREP.${NAME} \
                --time=10-00:00:00 \
                --output=${LOG_DIR}/_coverage_report_${NAME}.stdout \
                --error=${LOG_DIR}/_coverage_report_${NAME}.stderr \
                --wrap=\"bash ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/generate_coverage_report_single.sh \
                  $BAM
                  $COVERAGE_REPORT_FILE \
                  $MANIFEST\""

    echo "Current Dir: "$(pwd)
		echo $CMD
		if [ -f ${BUFFER_DIR}/nohup100.out ]; then
		  rm ${BUFFER_DIR}/nohup100.out
		fi
		echo $CMD >> ${BUFFER_DIR}/nohup100.out
		eval $CMD
	fi
done
