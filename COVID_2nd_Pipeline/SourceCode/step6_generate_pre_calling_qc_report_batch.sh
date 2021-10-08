#!/bin/sh

SCRIPT=$(readlink -f "$0")
DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc
#. ./global_config_bash.rc

MANIFEST_FILE=$1
dirBuildBAM=$2
strKTName=$3

DIRBAMRoot="/data/COVID_WGS/UpstreamAnalysis/PostPrimaryRun/Data/BAM/Batch"

#2 or 3 input (MANIFEST_FILE, dirBuildBAM and (PREQC_REPORT_FILE -- option))

DATE=`echo $(date +%m%d%Y)`
LOG_DIR=${CLUSTER_JOB_LOG_DIR}/${DATE}_${strKTName}

mkdir -p $LOG_DIR 2>/dev/null
if [ -f ${BUFFER_DIR}/nohup100.out ]; then
  rm ${BUFFER_DIR}/nohup100.out
fi

#create BAM list file
strBAMList=${dirBuildBAM}/bam.lst
if [ -f "${strBAMList}" ]; then
  rm ${strBAMList}
fi
# Create bam.list
find ${dirBuildBAM} -iname '*.bam' > ${strBAMList}

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

if [[ $# == 4 ]]; then
	# if there is an argument, enter patch mode
	PATCH_MODE=true
	PREQC_REPORT_FILE=$4
	echo "Entering patch mode to fill up the \"holes\" in existing coverage report file $PREQC_REPORT_FILE"
else
	# otherwise automatically generate the report file name
	PATCH_MODE=false
	PREQC_REPORT_FILE=${COVERAGE_REPORT_DIR}/pre_calling_qc_report_${DATE}_${strKTName}.txt
	
#	REPORT_LINE="Analysis ID\tAVERAGE_READ_LENGTH\tPERCENT_TOTAL_READS_WITH_3LESS_CUT\tPF_HQ_MEDIAN_MISMATCHES\tPF_HQ_MISMATCH_RATE_F\tPF_HQ_MISMATCH_RATE_R\tSTRAND_BALANCE_F\tSTRAND_BALANCE_R\tBAD_CYCLES\tPCT_CHIMERAS\tPCT_ADAPTER"\
#"\tPERCENT_PF_HQ_ALIGNED_READS\tPERCENT_PF_HQ_ALIGNED_BASES_L\tPERCENT_PF_HQ_ALIGNED_Q20_BASES_L\tPERCENT_PF_HQ_ALIGNED_BASES_R\tPERCENT_PF_HQ_ALIGNED_Q20_BASES_R\tBASES_QUAL_AVE"\
#"\tAT_DROPOUT\tGC_DROPOUT\tHS_AT_DROPOUT\tHS_GC_DROPOUT\tFOLD_80_BASE_PENALTY\tMEDIAN_INSERT_SIZES\tMEAN_INSERT_SIZES\tAVG_COVERAGE"\
#"\tLOWEST_PREADATPER_TOTAL_QSCORE_BASE\tLOWEST_PREADATPER_TOTAL_QSCORE\tLOWEST_BAITBIAS_TOTAL_QSCORE_BASE\tLOWEST_BAITBIAS_TOTAL_QSCORE\tOXIDATION_Q_AVE\tLOWEST_OXIDATION_Q\tLOWEST_OXIDATION_Q_CONTEXT"\
#"\tPERCENT_LOW_COVERAGE_CAPTURE_LOCI\tPERCENT_NO_COVERAGE_CAPTURE_LOCI\tPERCENT_POOR_MAPPING_QUALITY_LOCI\tPERCENT_READS_ENDING_IN_INDELS\tPERCENT_READS_ENDING_IN_SOFTCLIPS\tMEDIAN_INSERT_SIZES\tMEAN_INSERT_SIZES"

REPORT_LINE="Analysis ID\tAVERAGE_READ_LENGTH\tPERCENT_TOTAL_READS_WITH_3LESS_CUT\tPF_HQ_MEDIAN_MISMATCHES\tPF_HQ_MISMATCH_RATE_F\tPF_HQ_MISMATCH_RATE_R\tSTRAND_BALANCE_F\tSTRAND_BALANCE_R\tBAD_CYCLES\tPCT_CHIMERAS\tPCT_ADAPTER"\
"\tPERCENT_PF_HQ_ALIGNED_READS\tPERCENT_PF_HQ_ALIGNED_BASES_L\tPERCENT_PF_HQ_ALIGNED_Q20_BASES_L\tPERCENT_PF_HQ_ALIGNED_BASES_R\tPERCENT_PF_HQ_ALIGNED_Q20_BASES_R\tBASES_QUAL_AVE"\
"\tAT_DROPOUT\tGC_DROPOUT  \tMEAN_COVERAGE\tSD_COVERAGE\tMEDIAN_COVERAGE\tPCT_EXC_MAPQ\tPCT_EXC_UNPAIRED\tPCT_EXC_DUPE\tPCT_EXC_OVERLAP\tESTIMATED_LIBRARY_SIZE"\
"\tLOWEST_PREADATPER_TOTAL_QSCORE_BASE\tLOWEST_PREADATPER_TOTAL_QSCORE\tLOWEST_BAITBIAS_TOTAL_QSCORE_BASE\tLOWEST_BAITBIAS_TOTAL_QSCORE\tOXIDATION_Q_AVE"\
"\tPERCENT_LOW_COVERAGE_CAPTURE_LOCI\tPERCENT_NO_COVERAGE_CAPTURE_LOCI\tPERCENT_POOR_MAPPING_QUALITY_LOCI\tPERCENT_READS_ENDING_IN_INDELS\tPERCENT_READS_ENDING_IN_SOFTCLIPS\tMEDIAN_INSERT_SIZES\tMEAN_INSERT_SIZES"


#L="${NAME}\t${AVERAGE_READ_LENGTH}\t${PERCENT_TOTAL_READS_WITH_3LESS_CUT}\t${PF_HQ_MEDIAN_MISMATCHES}\t${PF_HQ_MISMATCH_RATE_F}\t${PF_HQ_MISMATCH_RATE_R}\t${STRAND_BALANCE_F}\t${STRAND_BALANCE_R}\t${BAD_CYCLES}\t${PCT_CHIMERAS}\t${PCT_ADAPTER}"\
#"\t${PERCENT_PF_HQ_ALIGNED_READS}\t${PERCENT_PF_HQ_ALIGNED_BASES_L}\t${PERCENT_PF_HQ_ALIGNED_Q20_BASES_L}\t${PERCENT_PF_HQ_ALIGNED_BASES_R}\t${PERCENT_PF_HQ_ALIGNED_Q20_BASES_R}\t${BASES_Q_AVE}"\
#"\t${AT_DROPOUT}\t${GC_DROPOUT}  \t${MEAN_COVERAGE}\t${SD_COVERAGE}\t${MEDIAN_COVERAGE}\t${PCT_EXC_MAPQ}\t${PCT_EXC_UNPAIRED}\t${PCT_EXC_DUPE}\t${PCT_EXC_OVERLAP}\t${ESTIMATED_LIBRARY_SIZE}"\
#"\t${LOWEST_PREADATPER_TOTAL_QSCORE_BASE}\t${LOWEST_PREADATPER_TOTAL_QSCORE}\t${LOWEST_BAITBIAS_TOTAL_QSCORE_BASE}\t${LOWEST_BAITBIAS_TOTAL_QSCORE}\t${OXIDATION_Q_AVE}"\
#"\t${PERCENT_LOW_COVERAGE_CAPTURE_REGION}\t${PERCENT_NO_COVERAGE_CAPTURE_REGION}\t${PERCENT_POOR_MAPPING_QUALITY}\t${PERCENT_READS_ENDING_IN_INDELS}\t${PERCENT_READS_ENDING_IN_SOFTCLIPS}\t${MEDIAN_INSERT_SIZES}\t${MEAN_INSERT_SIZES}"

  echo -e $REPORT_LINE > $PREQC_REPORT_FILE
fi	
	

# all new BAM files are generated under the incoming directory
# these are typically the files that need coverage report
echo "awk -F "\t" -v col='GROUP' '{for (i=1;i<=NF;i++) if ($i==col) print i; exit}' $MANIFEST_FILE"
GROUP_COL=$(awk -F "\t" -v col='GROUP' '{for (i=1;i<=NF;i++) if ($i==col) print i; exit}' $MANIFEST_FILE)
echo $GROUP_COL
ANALYSISID_COL=$(awk -F "\t" -v col='ANALYSIS ID' '{for (i=1;i<=NF;i++) if ($i~col) print i; exit}' $MANIFEST_FILE)
echo $ANALYSISID_COL
#BAM_DIR=$(dirname ${MANIFEST_FILE})/../bam_location
BAM_DIR=${dirBuildBAM}
cd $BAM_DIR

#create flag dir
strFlagDir=${DIRBAMRoot}/${strKTName}/Flag/PreCallingQCReport
if [ ! -d ${strFlagDir} ]; then
  mkdir -p ${strFlagDir}
fi

#for BAM in $(awk -F "\t" -v group=$GROUP_COL -v analysisid=$ANALYSISID_COL '{if(NR>1){print $group"/"$group"_"$analysisid".bam"}}' $MANIFEST_FILE | sort | uniq); do
#		echo $BAM
for BAM in `cat ${strBAMList}`; do
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
#		CMD="qsub -q $QUEUE -cwd \
#			-o ${LOG_DIR}/_pre_calling_qc_report_${NAME}.stdout -e ${LOG_DIR}/_pre_calling_qc_report_${NAME}.stderr \
#			-N PRECALLING_QC.$NAME \
#			-S /bin/sh ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/pre_calling_qc_single.sh $BAM $PREQC_REPORT_FILE $MANIFEST_FILE"
    strFlagWorking="${strFlagDir}/_pre_calling_QC_report_${NAME}.flag.working"
    strFlagDone="${strFlagDir}/_pre_calling_QC_report_${NAME}.flag.done"
    if [[ ! -f ${strFlagWorking} ]] && [[ ! -f ${strFlagDone} ]]; then
      # init working flag
      touch ${strFlagWorking}

      #strBashScript="${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/pre_calling_qc_single.sh"
      strBashScript="${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/pre_calling_wgsqc_single.sh"

      CMD="sbatch --ntasks=8 \
                  --nodes=1 \
                  --job-name=PRECALLING_QC.${NAME} \
                  --time=10-00:00:00 \
                  --mem=20G \
                  --output=${LOG_DIR}/_pre_calling_qc_report_${NAME}.stdout \
                  --error=${LOG_DIR}/_pre_calling_qc_report_${NAME}.stderr \
                  --wrap=\"bash ${strBashScript} \
                            $BAM \
                            $PREQC_REPORT_FILE \
                            $MANIFEST_FILE \
                            ${strFlagWorking} \
                            ${strFlagDone}\""

      echo $CMD
      echo $CMD >> ${BUFFER_DIR}/nohup100.out
      eval $CMD
    fi
	fi
done
