#!/bin/bash
set -o pipefail

echo "***************************"
echo "***** Slurm Job Info ******"
echo "***************************"
echo "Job id is   : ${SLURM_JOB_ID}"
echo "Job name is : ${SLURM_JOB_NAME}"
echo "Node list   : ${SLURM_JOB_NODELIST}"
echo "Partition   : ${SLURM_JOB_PARTITION}"
echo "QOS         : ${SLURM_JOB_QOS}"
echo "Cores       : ${SLURM_CPUS_ON_NODE}"
echo "Date        : $(date)"
echo "***************************"
echo

if [ -z $ILLUMINA_PROCESSED_DATA_ROOT_DIR ]; then
  echo "Error: environment variable ILLUMINA_PROCESSED_DATA_ROOT_DIR is unset."
  exit 1
fi

dirCurScript=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

. /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/global_config_bash.rc

#takes in arguments
argIndex=1
for arg in "$@"
do
  if [[ ${argIndex} == 1 ]]; then
    RUN_ID=${arg}
  elif [[ ${argIndex} == 2 ]]; then
    PROJECT_ID=${arg}
  elif [[ ${argIndex} == 3 ]]; then
    INPUT_BASE_ID=${arg}
  elif [[ ${argIndex} == 4 ]]; then
    OUTPUT_REPORT_FILE=${arg}
  elif [[ ${argIndex} == 5 ]]; then
    BASE_SAMPLE_DIR=${arg}
  elif [[ ${argIndex} == 6 ]]; then
    RUNINFO_LANES=${arg}
  elif [[ ${argIndex} == 7 ]]; then
    JOB_NAME=${arg}
  elif [[ ${argIndex} == 8 ]]; then
    alignerTool=${arg}
  elif [[ ${argIndex} == 9 ]]; then
    refVersion=${arg}
  elif [[ ${argIndex} == 10 ]]; then
    logDir=${arg}
  elif [[ ${argIndex} == 11 ]]; then
    bamDir=${arg}
  elif [[ ${argIndex} == 12 ]]; then
    reportDir=${arg}
  elif [[ ${argIndex} == 13 ]]; then
    sampleFlagDir=${arg}
  elif [[ ${argIndex} == 14 ]]; then
    containUMI=${arg}
  elif [[ ${argIndex} == 15 ]]; then
    sampleDir=${arg}
  elif [[ ${argIndex} == 16 ]]; then
    flagWorking=${arg}
  elif [[ ${argIndex} == 17 ]]; then
    flagDone=${arg}
  fi
  argIndex=$((argIndex + 1))
done


# splice input base id into three IDs
# the sample ID is the 1st one to the 3rd one from the last (this trick allows "_" in the sample ID) sample ID is just the SB number
SAMPLE_NAME=`echo $INPUT_BASE_ID | awk -F "_" 'sub(FS $NF,x)' | awk -F "_" 'sub(FS $NF,x)'`
echo "${SAMPLE_NAME}"

# the barcode is the 2nd one from the last
BARCODE=`echo $INPUT_BASE_ID | awk -F "_" '{print $(NF-1)}'`
echo "${BARCODE}"

# The sample name used in side CASAVA output: stats and report folder
sampleNameCASAVA="${SAMPLE_NAME}-${BARCODE}"
if [[ ${BARCODE} == "" ]]; then
  sampleNameCASAVA="${SAMPLE_NAME}"
fi
echo "${sampleNameCASAVA}"

# the lane id is the last one
LANE_ID=`echo $INPUT_BASE_ID | awk -F "_" '{print $NF}' | cut -c4-`

WORK_BASE_DIR=$ILLUMINA_PROCESSED_DATA_ROOT_DIR/$RUN_ID

#For getting log from trimming  -> Go after lunch
logMainDir=${WORK_BASE_DIR}/${ILLUMINA_LOG_SUB_DIR}


echo ====================
date
echo Report Generation
echo Run ID: $RUN_ID
echo Project ID: $PROJECT_ID
echo Sample Token: $INPUT_BASE_ID
echo Append Stats To: $OUTPUT_REPORT_FILE
echo ====================


##################################
# Report on sequencing quality
##################################
# first search for the logs for quality trimming
# and extract quality trimming stats
#LOG_FILE=${logMainDir}/*_quality_trimming_${INPUT_BASE_ID}.stderr
#echo "Extracting sequencing statistics from log file: $LOG_FILE"
#STAT_LINE=`grep "Input Read Pairs" $LOG_FILE | head -1`
#READ_PAIRS=`echo $STAT_LINE | awk '{print $4}'`
#let READS=$READ_PAIRS*2
#HQ_PAIRS=`echo $STAT_LINE | awk '{print $7}'`
#let HQ_READS=$HQ_PAIRS*2
#echo "Total Reads: $READS; Total HQ Reads: $HQ_READS"


##################################
# Report on duplicates
##################################
# now search in the BAM filtering report file for more stats
REPORT_FILE="${reportDir}/${INPUT_BASE_ID}_HQ_paired_flagstats.txt"
echo "Extracting reference mapping flagstats from file: ${REPORT_FILE}"
MAPPED_READS_ORIGINAL=$(grep 'in Original BAM' ${REPORT_FILE} | awk -F ':' '{print $2}')
MAPPED_READS_DEDUP=$(grep 'in Dedup BAM' ${REPORT_FILE} | awk -F ':' '{print $2}')
PERCENT_DUP=`echo "($MAPPED_READS_ORIGINAL - $MAPPED_READS_DEDUP) / $MAPPED_READS_ORIGINAL" | bc -l | xargs -I {} printf "%5.4f" {}`

echo "Reads Mapped: $MAPPED_READS_ORIGINAL; After Duplicate Removal: $MAPPED_READS_DEDUP; Percent Dup: $PERCENT_DUP"

# we also wanted to separate PCR dups from optical dups
# I am using Picard to get the number of optical dup here
# since the total number of dups are already calculated from above
# Note a little tricky thing here: Picard may not use 100% of reads in the BAM file for library estimation (it appears to remove some unreasonably aligned 
# reads), and hence in this report file, the READ_PAIRS_EXAMINED and READ_PAIR_DUPLICATES do not exactly match what we got from above
# they are however very close. Hence we'll assume the READ_PAIR_OPTICAL_DUPLICATES report is the accurate "estimation" here.

#echo "Estimating library complexity..."
#BAM_ORIGINAL=${WORK_BASE_DIR}/${ILLUMINA_BAM_SUB_DIR}/${INPUT_BASE_ID}_${SUFFIX_HQ_PAIRED}_${SUFFIX_BAM_ORIGINAL}.bam
#REPORT_FILE=${WORK_BASE_DIR}/${ILLUMINA_REPORT_SUB_DIR}/${INPUT_BASE_ID}_${SUFFIX_BAM_ORIGINAL}_library_complexity.txt
# run picard
##### $PICARD_ESTIMATE_LIBRARY_COMPLEXITY INPUT=$BAM_ORIGINAL OUTPUT=$REPORT_FILE

#if [ ! -f $REPORT_FILE ]; then
#  echo "Output report from Picard EstimateLibraryComplexity $REPORT_FILE expected but not found. Broken Picard pipe?" 1>&2
#  exit 1
#fi
#chmod g+rw $REPORT_FILE

# JH Edit: looks like during initial reference mapping stage, picard markduplicates has regenerated duplication metric report 
# which shows slightly different dup estimation but almost the same
# we are going to use that file instead (to save time cost)

REPORT_FILE=${reportDir}/${INPUT_BASE_ID}_HQ_paired_duplication_metrics.txt
# this number is reported in the 8th line, 7th column of the output file
READ_PAIRS_OPTICAL_DUP=$(awk -F '\t' '{ if (NR==8) print $8 }' $REPORT_FILE)

# convert pairs to reads
let READS_OPTICAL_DUP=2*$READ_PAIRS_OPTICAL_DUP
PERCENT_OPTICAL_DUP=`echo "$READS_OPTICAL_DUP / $MAPPED_READS_ORIGINAL" | bc -l | xargs printf "%5.4f" `

numPairedDup=$(awk -F '\t' '{ if (NR==8) print $7 }' ${REPORT_FILE})
let numTotalDup=2*${numPairedDup}
let READS_PCR_DUP=${numTotalDup}-${READS_OPTICAL_DUP}

PERCENT_PCR_DUP=`echo "$READS_PCR_DUP / $MAPPED_READS_ORIGINAL" | bc -l | xargs printf "%5.4f" `
echo "Reads PCR Dup: $READS_PCR_DUP ($PERCENT_PCR_DUP); Reads Optical Dup: $READS_OPTICAL_DUP ($PERCENT_OPTICAL_DUP) "

##################################
# Report on insertsize
##################################

REPORT_FILE=${reportDir}/${INPUT_BASE_ID}_HQ_paired_insertsize_metrics.txt

#get the median and mean insert sizes
MEDIAN_INSERT_SIZES=`awk -F '\t' '{ if (NR==8) print $1 }' $REPORT_FILE`
MEAN_INSERT_SIZES=`awk -F '\t' '{ if (NR==8) print $6 }' $REPORT_FILE`

##################################
# Report on coverage
##################################
# generate the coverage stats for the high-quality, dedupped, properly paired
BAM_FILE=${bamDir}/${INPUT_BASE_ID}_HQ_paired_dedup_nophix.bam

#### coverage report based on Illumina capture kit ####

REPORT_FILE=${reportDir}/${INPUT_BASE_ID}_HQ_paired_dedup_nophix_coverage_illumina.txt
#BED_FILE=$ILLUMINA_CAPTURE_KIT_BED_FILE
 
# generate the base-by-base pileup based on the designated bed file, extract the most essential 
# information (chromosome, position, read depth) into the report file
echo "Building pileup for coverage analysis based on Illumina Capture Kit..."
####

if [[ ! -f $REPORT_FILE ]]; then
#	samtools mpileup -l $BED_FILE $BAM_FILE | awk '{print $1"\t"$2"\t"$4}' > $REPORT_FILE
	samtools mpileup $BAM_FILE | awk '{print $1"\t"$2"\t"$4}' > $REPORT_FILE
fi
####

if [[ ! -f $REPORT_FILE ]]; then
  echo "Output report from samtools mpileup $REPORT_FILE expected but not found." 1>&2
  exit 1
fi
chmod g+rw $REPORT_FILE

totalNumBaseV37=3101804739
totalNumBaseV38=3217346917

TOTAL_BASES_COVERED=$(wc $REPORT_FILE -l | awk '{print $1}')
PERCENT_BASES_COVERED=$(echo $TOTAL_BASES_COVERED/${totalNumBaseV38} | bc -l  | xargs -I {} printf "%5.4f" {})


AVG_COVERAGE=$(awk -v total=${totalNumBaseV38} '{sum+=$3} END {print sum/total}' $REPORT_FILE | xargs printf "%2.2f")
echo "Reference to Whole genome, Total Exome Bases Covered: ${TOTAL_BASES_COVERED}; Percent: ${PERCENT_BASES_COVERED}; Average Coverage: ${AVG_COVERAGE}"

# CMD="awk -v threshold=5 '{if(\$3>=threshold) L++} END {print L}' $REPORT_FILE"
# TOTAL_BASES_5X=`eval $CMD`
# PERCENT_BASES_5X=`echo $TOTAL_BASES_5X/3101804739 | bc -l  | xargs -I {} printf "%5.4f" {}` 
# echo "Reference to Whole genome, Total Exome Bases >= 5x Coverage: ${TOTAL_BASES_5X}; Percent: ${PERCENT_BASES_5X}"

CMD="awk -v threshold=10 '{if(\$3>=threshold) L++} END {print L}' $REPORT_FILE"
TOTAL_BASES_10X=`eval $CMD`
PERCENT_BASES_10X=`echo $TOTAL_BASES_10X/${totalNumBaseV38} | bc -l  | xargs -I {} printf "%5.4f" {}`
echo "Reference to Whole genome, Total Exome Bases >= 10x Coverage: ${TOTAL_BASES_10X}; Percent: ${PERCENT_BASES_10X}"

CMD="awk -v threshold=15 '{if(\$3>=threshold) L++} END {print L}' $REPORT_FILE"
TOTAL_BASES_15X=`eval $CMD`
PERCENT_BASES_15X=`echo $TOTAL_BASES_15X/${totalNumBaseV38} | bc -l  | xargs -I {} printf "%5.4f" {}`
echo "Reference to Whole genome, Total Exome Bases >= 15x Coverage: ${TOTAL_BASES_15X}; Percent: ${PERCENT_BASES_15X}"

##Check UMI type -->
#casavaLaneDir=${WORK_BASE_DIR}/${ILLUMINA_CASAVA_SUB_DIR}/L${LANE_ID}
#TRIM_STATS=${casavaLaneDir}/Stats/AdapterTrimming.txt
#laneAllXML=${casavaLaneDir}/Reports/html/*/default/all/all/lane.html
## -->
#laneSampleXML=${casavaLaneDir}/Reports/html/*/*/${sampleNameCASAVA}/all/lane.html
## <--
#
#read1Caption="Read: 1"
#read1Index=1
#read2Caption="Read: 2"
#read2Index=2
#if [[ ${containUMI} == "yes" ]]; then
#  read2Caption="Read: 3"
#  read2Index=3
#fi
##<--
#echo
#echo "==== TRIM STATS File ===="
#echo ${TRIM_STATS}
#echo
#
#FORWARD_TRIMMED_BASE=`awk -F "\t" -v SAMPLE=${sampleNameCASAVA} -v ReadsIndex=${read1Index} -v RUNINFO_LANES=$RUNINFO_LANES '{ if ($2==ReadsIndex && $5==SAMPLE) {sum+=$8}} END {print sum/RUNINFO_LANES}' $TRIM_STATS`
#BACKWARD_TRIMMED_BASE=`awk -F "\t" -v SAMPLE=${sampleNameCASAVA} -v ReadsIndex=${read2Index} -v RUNINFO_LANES=$RUNINFO_LANES '{ if ($2==ReadsIndex && $5==SAMPLE) {sum+=$8}} END {print sum/RUNINFO_LANES}' $TRIM_STATS`
#
#FORWARD_TRIMMED_SEQ=`grep -A 36 "${read1Caption}" $TRIM_STATS | grep -A 36 -w "${sampleNameCASAVA}" | grep -v "Lane:"|awk -F "\t" -v RUNINFO_LANES=$RUNINFO_LANES '{sum+=$2} END {print sum/RUNINFO_LANES}'`
#BACKWARD_TRIMMED_SEQ=`grep -A 36 "${read2Caption}" $TRIM_STATS | grep -A 36 -w "${sampleNameCASAVA}" | grep -v "Lane:"|awk -F "\t" -v RUNINFO_LANES=$RUNINFO_LANES '{sum+=$2} END {print sum/RUNINFO_LANES}'`
#
#PERCENT_BIN_UNINDICES=0
#if [ -f ${laneAllXML} ]; then
#  PERCENT_BIN_UNINDICES=$(head -n 41 ${laneAllXML} | tail -n 1 | cut -d\< -f2 |  cut -d\> -f2)
#fi
#PERCENT_OF_LANE=0
#if [ -f ${laneSampleXML} ]; then
#  PERCENT_OF_LANE=$(head -n 41 ${laneSampleXML} | tail -n 1 | cut -d\< -f2 |  cut -d\> -f2)
#fi
#PERCENT_BIN_INDICES=$(echo 100-$PERCENT_BIN_UNINDICES | bc -l | xargs printf "%5.4f")

# Calculate Dimer Ratio -->
# (1) prepare folder/file path
resultFileName=L${LANE_ID}_${BASE_SAMPLE_DIR}_dimer_ratio
resultFilePath=${reportDir}/$resultFileName
# (2) Call python function of CalcDimerRatio
CMD="python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/CalAdapterDimerSingle.py ${sampleDir} $resultFilePath"
echo $CMD
eval $CMD
#Check return value and determine if we need to exit the main code
if [[ $? -ne 0 ]] || [[ ! -f $resultFilePath ]]; then
  echo "Error: CalAdapterDimerSingle.py was failed"
  exit 1
fi
# (3) Get calculation result
readsNum=$(grep -A 1 "Result:" ${resultFilePath} | tail -n 1 | awk -F '\t' '{print $1}' | sed 's/^[ \t]*//g')
dimerNum=$(grep -A 1 "Result:" ${resultFilePath} | tail -n 1 | awk -F '\t' '{print $2}' | sed 's/^[ \t]*//g')
dimerRatio=$(grep -A 1 "Result:" ${resultFilePath} | tail -n 1 | awk -F '\t' '{print $3}' | sed 's/^[ \t]*//g' | awk -F '%' '{print $1}')
#<----

# now print out all stats to the report file
# make sure the header of this file reflects the same order
# the header is being printed by another script generate_report_batch.sh

L="${SAMPLE_NAME},${BARCODE},${LANE_ID},${PROJECT_ID},${READS},"\
"${readsNum},${MAPPED_READS_ORIGINAL},${MAPPED_READS_DEDUP},${READS_PCR_DUP},${PERCENT_PCR_DUP},${READS_OPTICAL_DUP},${PERCENT_OPTICAL_DUP},${PERCENT_DUP},"\
"${MAPPED_READS_DEDUP},,,,${MEDIAN_INSERT_SIZES},${MEAN_INSERT_SIZES},"\
"${TOTAL_BASES_COVERED},${PERCENT_BASES_COVERED},${AVG_COVERAGE},${TOTAL_BASES_10X},${PERCENT_BASES_10X},${TOTAL_BASES_15X},${PERCENT_BASES_15X},,,,,,,,,"\
"${FORWARD_TRIMMED_BASE},${BACKWARD_TRIMMED_BASE},${FORWARD_TRIMMED_SEQ},${BACKWARD_TRIMMED_SEQ},${PERCENT_OF_LANE},${PERCENT_BIN_INDICES},"\
"${dimerNum},${dimerRatio}"

echo "Contents to be streamed to the flowcell-level report file:"
echo -e $L
#echo -e $L >> $OUTPUT_REPORT_FILE

# set flags
echo "Setting done flag..."
# set flags
if [[ -f ${flagWorking} ]];  then
  rm ${flagWorking}
fi
touch ${flagDone}

echo ====================
date
echo Done!
echo ====================
