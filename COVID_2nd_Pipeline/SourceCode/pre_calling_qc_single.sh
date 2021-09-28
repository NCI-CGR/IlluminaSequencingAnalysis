#!/bin/sh
# source global configuration files
# SCRIPT=$(readlink -f "$0")
# DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
# echo ${DCEG_SEQ_POOL_SCRIPT_DIR:-.} #this one doesn't work and print out to be this path /cm/local/apps/sge/var/spool/node115/job_scripts
# . ${DCEG_SEQ_POOL_SCRIPT_DIR}/global_config_bash.rc
. ./global_config_bash.rc
module load R/3.4.3
module load jdk/1.8.0_111
which java
IN_BAM=$1
OUT_REPORT=$2
MANIFEST=$3

OUT_DIR=${BUFFER_DIR}/PRE_QC
CUSTOMIZED_BED_DIR=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit
REF_DICT=$(echo $REFERENCE_GENOME | cut -f1 -d.).dict

NAME=$(basename $IN_BAM .bam)
#JAVA=/DCEG/Resources/Tools/jdk/7.17/jdk1.7.0_17/bin/java
#GATK=/DCEG/Projects/Exome/SequencingData/GATK_binaries/Latest/GenomeAnalysisTK.jar
#PICARD=/DCEG/Resources/Tools/Picard/Picard-2.10.10/picard.jar

#rm -rf ${BUFFER_DIR}/PRE_QC/${BUFFER_DIR}/
mkdir -p ${BUFFER_DIR}/PRE_QC/${NAME}/ 2>/dev/null

#retrieve captureKit from manifest file and use the corresponding bed file
ANALYSIS_ID=$(echo $NAME |cut -f2- -d_ | sed 's/_bqsr_final_out//g')
echo $ANALYSIS_ID
ASSAYID_COL=$(awk -F "\t" -v col='ASSAYID' '{for (i=1;i<=NF;i++) if ($i==col) print i; exit}' $MANIFEST)
#CAPTUREKIT=`awk -F '\t' -v QUERY_SAMPLE=$LIMSID '{if ($8 == QUERY_SAMPLE) {print $10}}' $MANIFEST`
CAPTUREKIT=`grep $ANALYSIS_ID $MANIFEST | awk -F '\t' -v col=$ASSAYID_COL '{print $col}' `
echo $MANIFEST
CAPTUREKIT=`echo $CAPTUREKIT | cut -d" " -f1`
echo $CAPTUREKIT
if [[ $CAPTUREKIT == *3.0* ]]; then
	PICARD_BAIT_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/SeqCap_EZ_Exome_v3_capture.interval_list
	PICARD_TARGET_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/SeqCap_EZ_Exome_v3_primary.interval_list
	GATK_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/SeqCap_EZ_Exome_v3_capture.intervals
	TOTAL_BASE=$EXOME_TARGETS_TOTAL_BASES_v3
	EXOME_TARGETS=$EXOME_TARGETS_v3
elif [[ $CAPTUREKIT == *UTR* ]]; then
    PICARD_BAIT_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/120430_HG19_ExomeV3_UTR_EZ_HX1_capture.interval_list
	PICARD_TARGET_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/120430_HG19_ExomeV3_UTR_EZ_HX1_primary.interval_list
	GATK_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/120430_HG19_ExomeV3_UTR_EZ_HX1_capture.intervals
	TOTAL_BASE=$EXOME_TARGETS_TOTAL_BASES_v3plusUTR
	EXOME_TARGETS=$EXOME_TARGETS_v3plusUTR
elif [[ $CAPTUREKIT == *2.0* ]]; then
	PICARD_BAIT_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/SeqCap_EZ_Exome_v3_capture.interval_list
	PICARD_TARGET_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/SeqCap_EZ_Exome_v3_primary.interval_list
	GATK_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/SeqCap_EZ_Exome_v3_capture.intervals
	TOTAL_BASE=$EXOME_TARGETS_TOTAL_BASES_v3
	EXOME_TARGETS=$EXOME_TARGETS_v3
elif [[ $CAPTUREKIT == *EZ_Choice_Kid-Lung-Extra* ]]; then
    PICARD_BAIT_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/161006_HG19_Kid-Lung-Extra_EZ_HX3_capture_targets.interval_list
	PICARD_TARGET_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/161006_HG19_Kid-Lung-Extra_EZ_HX3_primary_targets.interval_list
	GATK_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/161006_HG19_Kid-Lung-Extra_EZ_HX3_primary_targets.intervals
	TOTAL_BASE=1780492
	EXOME_TARGETS=${CUSTOMIZED_BED_DIR}/161006_HG19_Kid-Lung-Extra_EZ_HX3_primary_targets.bed
elif [[ $CAPTUREKIT == *Karina-XP* ]]; then
        PICARD_BAIT_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/Karina-XP_capture_targets.interval_list
	PICARD_TARGET_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/Karina-XP_primary_targets.interval_list
	GATK_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/Karina-XP_primary_targets.intervals
	TOTAL_BASE=941427
	EXOME_TARGETS=${CUSTOMIZED_BED_DIR}/Karina-XP_primary_targets.bed	
elif [[ $CAPTUREKIT == *Dean-Koshiol-4* ]]; then
    PICARD_BAIT_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/Dean-Koshiol-4_capture.interval_list
	PICARD_TARGET_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/Dean-Koshiol-4_primary.interval_list
	GATK_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/Dean-Koshiol-4_capture.intervals
	TOTAL_BASE=2956332
	EXOME_TARGETS=/DCEG/CGF/Laboratory/LIMS/Oligo_Orders/NimblegenCapture/Dean-Koshiol-4/Selection_Results/Dean-Koshiol-4_capture_targets.bed
elif [[ $CAPTUREKIT == *Agilent_exome_test* ]]; then
    PICARD_BAIT_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/S31285117_Regions.interval_list
	PICARD_TARGET_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/S31285117_Regions.interval_list
	GATK_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/S31285117_Regions.intervals
	TOTAL_BASE=35804808
	EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/S31285117_Regions.bed
elif [[ $CAPTUREKIT == *Roche_exome_test* ]]; then
    PICARD_BAIT_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/New_RSSExome_capture_targets.hg19.interval_list
	PICARD_TARGET_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/New_RSSExome_primary_targets.hg19.interval_list
	GATK_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/New_RSSExome_capture_targets.hg19.intervals
	TOTAL_BASE=42836424
	EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/New_RSSExome_capture_targets.hg19.bed
elif [[ $CAPTUREKIT == *IDT_exome_test* ]]; then
    PICARD_BAIT_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/xgen-exome-research-panel-targets.interval_list
	PICARD_TARGET_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/xgen-exome-research-panel-targets.interval_list
	GATK_INTERVALS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/xgen-exome-research-panel-targets.intervals
	TOTAL_BASE=39086460
	EXOME_TARGETS=/DCEG/Projects/Exome/SequencingData/BED_FILES/customized_capturekit/xgen-exome-research-panel-targets.bed
else
	EXOME_TARGETS=$(ls ${CUSTOMIZED_BED_DIR}/*${CAPTUREKIT}*_capture.bed)
	PRIMARY_BED=$(ls ${CUSTOMIZED_BED_DIR}/*${CAPTUREKIT}*_primary.bed)

	if [[ ! -f ${EXOME_TARGETS} ]] || [[ ! -f ${PRIMARY_BED} ]]; then
	    echo no bed file found!
	    exit 0
	fi
	TOTAL_BASE=$( awk '{print $3-$2}' ${EXOME_TARGETS} | awk '{sum+=$1} END {print sum}')
    PICARD_BAIT_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/${CAPTUREKIT}_capture.interval_list
	PICARD_TARGET_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/${CAPTUREKIT}_primary.interval_list
	GATK_INTERVALS=${BUFFER_DIR}/INTERVAL_FILES/${CAPTUREKIT}_capture.intervals
	if [[ ! -f ${PICARD_BAIT_INTERVALS} ]]; then
	    CMD="java -jar $PICARD BedToIntervalList I=${EXOME_TARGETS} O=${BUFFER_DIR}/INTERVAL_FILES/${CAPTUREKIT}_capture.interval_list SD=${REF_DICT}"
	    echo $CMD
	    eval $CMD
	    awk -F"\t" '{print $1":"$2+1"-"$3}' ${EXOME_TARGETS} > ${GATK_INTERVALS}
    fi
    if [[ ! -f ${PICARD_TARGET_INTERVALS} ]]; then
	    
	    CMD="java -jar $PICARD BedToIntervalList I=${PRIMARY_BED} O=${BUFFER_DIR}/INTERVAL_FILES/${CAPTUREKIT}_primary.interval_list SD=${REF_DICT}"
	    echo $CMD
	    eval $CMD
    fi
fi
echo $EXOME_TARGETS

echo ====================

#if IN_BAM is newer than any of the qc metrics, then IN_BAM is newly created, redo the pre qc metrics generation
if [[ -f ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.gc_bias.summary_metrics && `find $IN_BAM  -printf "%l\n" | xargs -I {} find {} -newer ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.gc_bias.summary_metrics` ]]; then 
	echo ${IN_BAM} file newer than last created qc metrics, replace old metrics.
	REPLACE_TABLE=true
else 
	echo ${IN_BAM} already has most updated qc metrics, grep the stats directly.
	REPLACE_TABLE=false
fi

echo $REPLACE_TABLE

#picard CollectMultipleMetrics
#collect multiple classes of metrics. This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time to cut down on the time spent reading in data from input files. Available modules include CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, and CollectQualityYieldMetrics. The tool produces outputs of '.pdf' and '.txt' files for each module, except for the CollectAlignmentSummaryMetrics module, which outputs only a '.txt' file. Output files are named by specifying a base name (without any file extensions).


date
echo Running all tools to collect qc metrics.  
echo ====================

CMD="java -Xmx16g -jar $PICARD CollectMultipleMetrics I=${TMP_DIR}/${NAME}.bam O=${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics R=$REFERENCE_GENOME PROGRAM=null PROGRAM=CollectSequencingArtifactMetrics PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectGcBiasMetrics VALIDATION_STRINGENCY=LENIENT"
#CollectSequencingArtifactMetrics substitution_rate of the 6 base changes ,pre_adapter_summary_metrics TOTAL_QSCORE <35
#CollectAlignmentSummaryMetrics %PF_HQ_ALIGNED_READS(mapping quality of >20) %PF_HQ_ALIGNED_Q20_BASES PCT_READS_ALIGNED_IN_PAIRS PCT_PF_READS_IMPROPER_PAIRS PF_MISMATCH_RATE STRAND_BALANCE
#PF_NOISE_READS: number of PF reads that are composed entirly out of As and/or Ns
#PF_HQ_MEDIAN_MISMATCHES: median number of mismatches in PF_HQ_ALIGNED_READS
#PF_HQ_ERROR_RATE: percentage of bases that mismatch in PF_HQ_ALIGNED_READS
#BAD_CYCLES: Number of instrument cycles where 80% or more base calls were no-calls
#PCT_CHIMERAS: percentage of reads that map outside of a maximum insert size (usually 100kb) or that have the two ends mapping to different chromosomes
#PCT_ADAPTER: percentage of unaligned PF reads that match to a known adapter sequence right from the start of the read.
#CollectGcBiasMetrics AT_DROPOUT GC_DROPOUT

if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.gc_bias.summary_metrics ]] || [[ $REPLACE_TABLE != "false" ]]; then
	echo "java -Xmx16g -jar $PICARD ReorderSam I=$IN_BAM O=${TMP_DIR}/${NAME}.bam reference=$REFERENCE_GENOME CREATE_INDEX=true"
	java -Xmx16g -jar $PICARD ReorderSam I=$IN_BAM O=${TMP_DIR}/${NAME}.bam reference=$REFERENCE_GENOME CREATE_INDEX=true
	echo $CMD
	eval $CMD
fi

CMD="java -Xmx16g -jar /DCEG/Resources/Tools/Picard/Picard-2.10.10/picard.jar QualityScoreDistribution I=$IN_BAM O=${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt CHART=${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.pdf VALIDATION_STRINGENCY=LENIENT"
#QualityScoreDistribution can get the percentage of bases higher than 30

if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt ]] || [[ $REPLACE_TABLE != "false" ]]; then
	echo $CMD
	eval $CMD
fi

#Picard markduplicates has to be run and duplicate reads have to be retained before CollectHsMetrics to enable the HS_LIBRARY_SIZE field 
#CMD="java -Xmx16g -jar /DCEG/Resources/Tools/Picard/Picard-2.10.10/picard.jar MarkDuplicates I=$IN_BAM O=${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.bam M=${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.txt"

#if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.txt || $REPLACE_TABLE != "false" ]]; then
#	echo $CMD
#	eval $CMD
#fi

#Collects hybrid-selection (HS) metrics for a SAM or BAM file. This tool takes a SAM/BAM file input and collects metrics that are specific for sequence datasets generated through hybrid-selection. Hybrid-selection (HS) is the most commonly used technique to capture exon-specific sequences for targeted sequencing experiments such as exome sequencing; for more information, please see the corresponding GATK Dictionary entry. 
#not very useful, similar function with coverage report
#AT_DROPOUT GC_DROPOUT
#FOLD_80_BASE_PENALTY: measure of non-uniformity. Number of units of sequencing necessary to raise 80% of the bases to the mean coverage (<5, between 3 and 4 for a good exome set)
#HS_LIBRARY_SIZE empty...
CMD="java -Xmx16g -jar $PICARD CollectHsMetrics I=$IN_BAM O=${BUFFER_DIR}/PRE_QC/${NAME}/hs_metrics.txt BAIT_INTERVALS=$PICARD_BAIT_INTERVALS TARGET_INTERVALS=$PICARD_TARGET_INTERVALS R=$REFERENCE_GENOME VALIDATION_STRINGENCY=LENIENT"

if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/hs_metrics.txt ]] || [[ $REPLACE_TABLE != "false" ]]; then
	echo $CMD
	eval $CMD
fi

#CollectOxoGMetrics #average OXIDATION_Q

CMD="java -Xmx16g -jar $PICARD CollectOxoGMetrics I=$IN_BAM O=${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt R=$REFERENCE_GENOME VALIDATION_STRINGENCY=LENIENT"

if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt ]] || [[ $REPLACE_TABLE != "false" ]]; then
	echo $CMD
	eval $CMD
fi

CMD="java -Xmx16g -jar $GATK -T CallableLoci -R $REFERENCE_GENOME -I $IN_BAM -summary ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt -o ${BUFFER_DIR}/PRE_QC/${NAME}/callable_status.bed -L $GATK_INTERVALS"
#table.txt POOR_MAPPING_QUALITY (<10),LOW_COVERAGE(<4),NO_COVERAGE,still the callable_status.bed file can serve as a final vcf filter
if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt ]] || [[ $REPLACE_TABLE != "false" ]]; then
	echo $CMD
	eval $CMD
fi

# CMD="java -Xmx16g -jar $GATK -T ErrorRatePerCycle  -R $REFERENCE_GENOME -I $IN_BAM -o ${BUFFER_DIR}/PRE_QC/${NAME}/error_rates.gatkreport.txt"
# #not very useful, need to run across more samples to see the trend

# if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/error_rates.gatkreport.txt || $REPLACE_TABLE != "false" ]]; then
	# echo $CMD
	# eval $CMD
# fi

CMD="java -Xmx16g -jar $PICARD CollectInsertSizeMetrics \
  INPUT=$IN_BAM \
  OUTPUT=${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt \
  HISTOGRAM_FILE=${BUFFER_DIR}/PRE_QC/${NAME}/InsertSizeHist.pdf \
  DEVIATIONS=10.0 MINIMUM_PCT=0.05 \
  VALIDATION_STRINGENCY=LENIENT"
if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt ]] || [[ $REPLACE_TABLE != "false" ]];then
	echo $CMD
	eval $CMD
fi


CMD="java -Xmx16g -jar $GATK -T CountTerminusEvent -R $REFERENCE_GENOME -I ${TMP_DIR}/${NAME}_filtered.bam -o ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt -L $GATK_INTERVALS"
#can try to add the two stats to the report

if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt ]] || [[ $REPLACE_TABLE != "false" ]]; then
	#filter in bam proper aligned reads first
	samtools view -f 3 -b $IN_BAM > ${TMP_DIR}/${NAME}_filtered.bam
	samtools index ${TMP_DIR}/${NAME}_filtered.bam
	echo $CMD
	eval $CMD
	rm -rf ${TMP_DIR}/${NAME}_filtered.bam* 2> /dev/null 

fi


CMD="java -Xmx16g -jar $GATK -T ReadLengthDistribution  -R $REFERENCE_GENOME -I $IN_BAM -o ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl"

if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl ]] || [[ $REPLACE_TABLE != "false" ]]; then
	echo $CMD
	eval $CMD
fi



CMD="samtools depth -d 0 -b  $EXOME_TARGETS $IN_BAM  | awk '{print \$3}' > ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}.coverage"

if [[ ! -s ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}.coverage ]]; then
	echo $CMD
	samtools depth -d 0 -b $EXOME_TARGETS $IN_BAM  | awk '{print $3}' > ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}.coverage
fi

echo ====================
date
echo Start collecting stats. 
echo ====================

#get the lowest total_qscore of each base change from pre_adapter_summary_metrics and bait_bias_summary_metrics. 35 can be a warning threshold
#remove the trailing empty lines and find the lowest score from line 8
LOWEST_PREADATPER_TOTAL_QSCORE=$(awk '/^$/ {nlstack=nlstack "\n";next;} {printf "%s",nlstack; nlstack=""; print;}' ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.pre_adapter_summary_metrics | awk 'NR == 8 {min=$5} $5 <= min {min = $5} END { print min }')
LOWEST_PREADATPER_TOTAL_QSCORE_BASE=`awk -v score=$LOWEST_PREADATPER_TOTAL_QSCORE '{if($5==score){print $3">"$4}}' ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.pre_adapter_summary_metrics`
LOWEST_BAITBIAS_TOTAL_QSCORE=$(awk '/^$/ {nlstack=nlstack "\n";next;} {printf "%s",nlstack; nlstack=""; print;}' ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.bait_bias_summary_metrics | awk 'NR == 8 {min=$5} $5 <= min {min = $5} END { print min }')
LOWEST_BAITBIAS_TOTAL_QSCORE_BASE=`awk -v score=$LOWEST_BAITBIAS_TOTAL_QSCORE '{if($5==score){print $3">"$4}}' ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.bait_bias_summary_metrics`

#get the %PF_HQ_ALIGNED_READS(mapping quality of >20) %PF_HQ_ALIGNED_Q20_BASES PCT_READS_ALIGNED_IN_PAIRS PCT_PF_READS_IMPROPER_PAIRS(add if we decided to remove proper align filter) PF_MISMATCH_RATE STRAND_BALANCE from CollectAlignmentSummaryMetrics
PF_NOISE_READS=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $5}'`
PF_HQ_MEDIAN_MISMATCHES=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $12}'`
#PF_HQ_ERROR_RATE=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $14}'`
BAD_CYCLES=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $21}'`


TOTAL_READS_ALIGNED=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $6}'`
PF_HQ_ALIGNED_READS=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $9}'`
if [[ ! $TOTAL_READS_ALIGNED -gt 0 ]]; then	
	echo "Error. Total reads aligned is less than 0, something is wrong."
else
	PERCENT_PF_HQ_ALIGNED_READS=`echo "$PF_HQ_ALIGNED_READS / $TOTAL_READS_ALIGNED" | bc -l | xargs -I {} printf "%5.4f" {}`
fi

PF_ALIGNED_BASES_L=`grep ^FIRST_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $8}'`
PF_HQ_ALIGNED_BASES_L=`grep ^FIRST_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $10}'`
PF_HQ_ALIGNED_Q20_BASES_L=`grep ^FIRST_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $11}'`

PF_ALIGNED_BASES_R=`grep ^SECOND_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $8}'`
PF_HQ_ALIGNED_BASES_R=`grep ^SECOND_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $10}'`
PF_HQ_ALIGNED_Q20_BASES_R=`grep ^SECOND_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $11}'`

if [[ ! $PF_ALIGNED_BASES_L -gt 0 ]]; then	
	echo "Error. PF_ALIGNED_BASES_L is less than 0, something is wrong."
else
	PERCENT_PF_HQ_ALIGNED_BASES_L=`echo "$PF_HQ_ALIGNED_BASES_L / $PF_ALIGNED_BASES_L" | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_PF_HQ_ALIGNED_Q20_BASES_L=`echo "$PF_HQ_ALIGNED_Q20_BASES_L / $PF_ALIGNED_BASES_L" | bc -l | xargs -I {} printf "%5.4f" {}`
fi
if [[ ! $PF_ALIGNED_BASES_R -gt 0 ]]; then	
	echo "Error. PF_ALIGNED_BASES_R is less than 0, something is wrong."
else	
	PERCENT_PF_HQ_ALIGNED_BASES_R=`echo "$PF_HQ_ALIGNED_BASES_R / $PF_ALIGNED_BASES_R" | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_PF_HQ_ALIGNED_Q20_BASES_R=`echo "$PF_HQ_ALIGNED_Q20_BASES_R / $PF_ALIGNED_BASES_R" | bc -l | xargs -I {} printf "%5.4f" {}`
fi

PF_HQ_MISMATCH_RATE_F=`grep ^FIRST_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $14}'`
PF_HQ_MISMATCH_RATE_R=`grep ^SECOND_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $14}'`

STRAND_BALANCE_F=`grep ^FIRST_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $22}'`
STRAND_BALANCE_R=`grep ^SECOND_OF_PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $22}'`
PCT_CHIMERAS=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $23}'`
PCT_ADAPTER=`grep ^PAIR ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics | awk -F "\t" '{print $24}'`

#get AT_DROPOUT GC_DROPOUT from CollectGcBiasMetrics
#AT_DROPOUT and GC_DROPOUT metrics, which indicate the percentage of misaligned reads that correlate with low (%-GC is < 50%) or high (%-GC is > 50%) GC content respectively
AT_DROPOUT=`awk -F "\t" '{if(NR==8){print $6}}' ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.gc_bias.summary_metrics`
GC_DROPOUT=`awk -F "\t" '{if(NR==8){print $7}}' ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.gc_bias.summary_metrics`

#QualityScoreDistribution can get the percentage of bases higher than 32 and lower than 12 of total bases 
TOTAL_BASES=`awk 'NR>8' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt | awk -F "\t" '{sum+=$2} END {print sum}'`
#ROW22=$(awk '{if($1==22){print NR}}' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt)
#ROW37=$(awk '{if($1==37){print NR}}' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt)
BASES_Q_AVE=`awk -F"\t" '{if(NR>8) {sum1+=$2;sum2+=$1*$2}} END{print sum2/sum1}' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt `
#BASES_Q_22DOWN=`awk -v row22=$ROW22 'NR>8&&NR<row22' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt | awk -F "\t" '{sum+=$2} END {print sum}'`
#echo $ROW22 $ROW37
#echo $BASES_Q_37UP $BASES_Q_22DOWN
#echo $TOTAL_BASES
if [[ ! $TOTAL_BASES -gt 0 ]]; then	
	echo "Error. TOTAL_BASES is less than 0, something is wrong."
#else
#	PERCENT_BASES_Q_37UP=`echo "$BASES_Q_37UP / $TOTAL_BASES" | bc -l | xargs -I {} printf "%5.4f" {}`
#	PERCENT_BASES_Q_22DOWN=`echo "$BASES_Q_22DOWN / $TOTAL_BASES" | bc -l | xargs -I {} printf "%5.4f" {}`
fi

#Collects hybrid-selection (HS) metrics AT_DROPOUT GC_DROPOUT FOLD_80_BASE_PENALTY HS_LIBRARY_SIZE
HS_AT_DROPOUT=`awk -F "\t" '{if(NR==8){print $51}}' ${BUFFER_DIR}/PRE_QC/${NAME}/hs_metrics.txt`
HS_GC_DROPOUT=`awk -F "\t" '{if(NR==8){print $52}}' ${BUFFER_DIR}/PRE_QC/${NAME}/hs_metrics.txt`
FOLD_80_BASE_PENALTY=`awk -F "\t" '{if(NR==8){print $35}}' ${BUFFER_DIR}/PRE_QC/${NAME}/hs_metrics.txt`
#HS_LIBRARY_SIZE=`awk -F "\t" '{if(NR==8){print $44}}' ${BUFFER_DIR}/PRE_QC/${NAME}/hs_metrics.txt`

MEDIAN_INSERT_SIZES=`awk -F '\t' '{ if (NR==8) print $1 }' ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt`
MEAN_INSERT_SIZES=`awk -F '\t' '{ if (NR==8) print $5 }' ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt`

#CollectOxoGMetrics average OXIDATION_Q. 4 Osteo problematic samples are showing ~20 of this score.

OXIDATION_Q_AVE=`awk 'NR>7&&NR<24' ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt |  awk '{ sum += $12 } END { print sum /16 }'`
LOWEST_OXIDATION_Q=$(awk '/^$/ {nlstack=nlstack "\n";next;} {printf "%s",nlstack; nlstack=""; print;}' ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt | awk 'NR == 8 {min=$12} $12 <= min {min = $12} END { print min }')
LOWEST_OXIDATION_Q_CONTEXT=`awk -v score=$LOWEST_OXIDATION_Q '{if($12==score){print $3}}' ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt`

#GATK CallableLoci  POOR_MAPPING_QUALITY (<10),LOW_COVERAGE(<4),NO_COVERAGE
LOW_COVERAGE_CAPTURE_REGION=`awk '{if (NR==5){print $2}}' ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt`
NO_COVERAGE_CAPTURE_REGION=`awk '{if (NR==4){print $2}}' ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt`
POOR_MAPPING_QUALITY=`awk '{if (NR==7){print $2}}' ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt`
#echo "$LOW_COVERAGE_CAPTURE_REGION / $TOTAL_BASE "
PERCENT_LOW_COVERAGE_CAPTURE_REGION=`echo "$LOW_COVERAGE_CAPTURE_REGION / $TOTAL_BASE " | bc -l | xargs -I {} printf "%5.4f" {}`
PERCENT_NO_COVERAGE_CAPTURE_REGION=`echo "$NO_COVERAGE_CAPTURE_REGION / $TOTAL_BASE " | bc -l | xargs -I {} printf "%5.4f" {}`
PERCENT_POOR_MAPPING_QUALITY=`echo "$POOR_MAPPING_QUALITY / $TOTAL_BASE " | bc -l | xargs -I {} printf "%5.4f" {}`

#GATK CountTerminusEvent, this test has to be down on proper paried filtered bam. or will get error message of SAM/BAM file SAMFileReader is malformed: read does not have any bases, it's all hard clips
READS_ENDING_IN_INDELS=`head -1 ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt | awk -F ":" '{print $2}'`
READS_ENDING_IN_SOFTCLIPS=`head -2 ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt  | tail -1 | awk -F ":" '{print $2}'`

#GATK ReadLengthDistribution
TOTAL_READS=`awk 'NR>4' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl | awk '{for(i=2;i<=NF;i++) sum+=$i} END {print sum}'`
#TOTAL_BASES=`awk 'NR>4' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl | awk '{sum+=$2*$1} END {print sum}'`
#for bam file with multiple lanes(multiple RG), there will be a column of reads number for each RG.
TOTAL_BASES=`awk 'NR>4' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl | awk 't=0;{for(i=2;i<=NF;i++) t+=$i; sum+=t*$1} END {print sum}'`
#echo $TOTAL_BASES
TOTAL_READS_WITH_3LESS_CUT=`awk '/^$/ {nlstack=nlstack "\n";next;} {printf "%s",nlstack; nlstack=""; print;}' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl | tail -4 | awk '{for(i=2;i<=NF;i++) sum+=$i} END {print sum}'`

if [[ ! $TOTAL_READS -gt 0 ]]; then	
	echo "Error. TOTAL_READS is less than 0, something is wrong."
else
	AVERAGE_READ_LENGTH=`echo "$TOTAL_BASES / $TOTAL_READS" | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_TOTAL_READS_WITH_3LESS_CUT=`echo "$TOTAL_READS_WITH_3LESS_CUT / $TOTAL_READS " | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_READS_ENDING_IN_INDELS=`echo "$READS_ENDING_IN_INDELS / $TOTAL_READS " | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_READS_ENDING_IN_SOFTCLIPS=`echo "$READS_ENDING_IN_SOFTCLIPS / $TOTAL_READS " | bc -l | xargs -I {} printf "%5.4f" {}`

fi

AVG_COVERAGE=`awk -v total=$TOTAL_BASE '{sum+=$1} END {print sum/total}' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}.coverage | xargs printf "%2.2f" `

L="${NAME}\t${AVERAGE_READ_LENGTH}\t${PERCENT_TOTAL_READS_WITH_3LESS_CUT}\t${PF_HQ_MEDIAN_MISMATCHES}\t${PF_HQ_MISMATCH_RATE_F}\t${PF_HQ_MISMATCH_RATE_R}\t${STRAND_BALANCE_F}\t${STRAND_BALANCE_R}\t${BAD_CYCLES}\t${PCT_CHIMERAS}\t${PCT_ADAPTER}"\
"\t${PERCENT_PF_HQ_ALIGNED_READS}\t${PERCENT_PF_HQ_ALIGNED_BASES_L}\t${PERCENT_PF_HQ_ALIGNED_Q20_BASES_L}\t${PERCENT_PF_HQ_ALIGNED_BASES_R}\t${PERCENT_PF_HQ_ALIGNED_Q20_BASES_R}\t${BASES_Q_AVE}"\
"\t${AT_DROPOUT}\t${GC_DROPOUT}\t${HS_AT_DROPOUT}\t${HS_GC_DROPOUT}\t${FOLD_80_BASE_PENALTY}\t${MEDIAN_INSERT_SIZES}\t${MEAN_INSERT_SIZES}\t${AVG_COVERAGE}"\
"\t${LOWEST_PREADATPER_TOTAL_QSCORE_BASE}\t${LOWEST_PREADATPER_TOTAL_QSCORE}\t${LOWEST_BAITBIAS_TOTAL_QSCORE_BASE}\t${LOWEST_BAITBIAS_TOTAL_QSCORE}\t${OXIDATION_Q_AVE}\t${LOWEST_OXIDATION_Q}\t${LOWEST_OXIDATION_Q_CONTEXT}"\
"\t${PERCENT_LOW_COVERAGE_CAPTURE_REGION}\t${PERCENT_NO_COVERAGE_CAPTURE_REGION}\t${PERCENT_POOR_MAPPING_QUALITY}\t${PERCENT_READS_ENDING_IN_INDELS}\t${PERCENT_READS_ENDING_IN_SOFTCLIPS}"

echo "Doing the last piece of housekeeping..."
CMD="rm ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}.coverage ${TMP_DIR}/${NAME}.bam ${TMP_DIR}/${NAME}.bai ${BUFFER_DIR}/PRE_QC/${NAME}/callable_status.bed ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}.coverage"
echo "Command: $CMD"
eval $CMD

echo "Contents to be streamed to the flowcell-level report file:"
echo -e $L
#echo -e $L >> $OUT_REPORT
echo ====================
date
echo Done!
echo ====================

