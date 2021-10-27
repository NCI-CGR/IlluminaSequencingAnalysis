#!/bin/sh
set -o pipefail
# source global configuration files
#. ${BUFFER_DIR:-.}/global_config_bash.rc
SCRIPT=$(readlink -f "$0")
DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc

module load R/4.1
module load java/1.8.0_211
which java

IN_BAM=$1
OUT_REPORT=$2
MANIFEST=$3
strFlagWorking=$4
strFlagDone=$5

OUT_DIR=${BUFFER_DIR}/PRE_QC
NAME=$(basename $IN_BAM .bam)

#JAVA=/DCEG/Resources/Tools/jdk/7.17/jdk1.7.0_17/bin/java
##GATK=/DCEG/Projects/Exome/SequencingData/GATK_binaries/Latest/GenomeAnalysisTK.jar
#PICARD=/DCEG/Resources/Tools/Picard/Picard-2.10.10/picard.jar
#REFERENCE_GENOME=/DCEG/Projects/PopulationExome/EAGLE/b37/human_g1k_v37_decoy.fasta

#rm -rf ${BUFFER_DIR}/PRE_QC/${BUFFER_DIR}/
mkdir -p ${BUFFER_DIR}/PRE_QC/${NAME}/ 2>/dev/null

TOTAL_BASE=3101804739
#retrieve captureKit from manifest file and use the corresponding bed file
ANALYSIS_ID=$(echo $NAME |cut -f2- -d_)
#ASSAYID_COL=$(awk -F "\t" -v col='ASSAYID' '{for (i=1;i<=NF;i++) if ($i==col) print i; exit}' $MANIFEST)
#CAPTUREKIT=`awk -F '\t' -v QUERY_SAMPLE=$LIMSID '{if ($8 == QUERY_SAMPLE) {print $10}}' $MANIFEST`
#CAPTUREKIT=`grep $ANALYSIS_ID $MANIFEST | awk -F '\t' -v col=$ASSAYID_COL '{print $col}' `
#echo $MANIFEST
#CAPTUREKIT=`echo $CAPTUREKIT | cut -d" " -f1`
#echo $CAPTUREKIT
#if [[ $CAPTUREKIT == *3.0* ]]; then
#	PICARD_BAIT_INTERVALS=/home/luow2/20170817_qc_test/SeqCap_EZ_Exome_v3_capture.interval_list
#	PICARD_TARGET_INTERVALS=/home/luow2/20170817_qc_test/SeqCap_EZ_Exome_v3_primary.interval_list
#	GATK_INTERVALS=/home/luow2/20170817_qc_test/SeqCap_EZ_Exome_v3_capture.intervals
#elif [[ $CAPTUREKIT == *UTR* ]]; then
#    PICARD_BAIT_INTERVALS=/home/luow2/20170817_qc_test/120430_HG19_ExomeV3_UTR_EZ_HX1_capture.interval_list
#	PICARD_TARGET_INTERVALS=/home/luow2/20170817_qc_test/120430_HG19_ExomeV3_UTR_EZ_HX1_primary.interval_list
#	GATK_INTERVALS=/home/luow2/20170817_qc_test/120430_HG19_ExomeV3_UTR_EZ_HX1_capture.intervals
#else
#	echo no bed file found!
#	exit 0
#    PICARD_BAIT_INTERVALS=/home/luow2/20170817_qc_test/120430_HG19_ExomeV3_UTR_EZ_HX1_capture.interval_list
#	PICARD_TARGET_INTERVALS=/home/luow2/20170817_qc_test/120430_HG19_ExomeV3_UTR_EZ_HX1_primary.interval_list
#	GATK_INTERVALS=/home/luow2/20170817_qc_test/120430_HG19_ExomeV3_UTR_EZ_HX1_capture.intervals
#fi

#picard CollectMultipleMetrics
#collect multiple classes of metrics. This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time to cut down on the time spent reading in data from input files. Available modules include CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, and CollectQualityYieldMetrics. The tool produces outputs of '.pdf' and '.txt' files for each module, except for the CollectAlignmentSummaryMetrics module, which outputs only a '.txt' file. Output files are named by specifying a base name (without any file extensions).

echo ====================
date
echo Running all tools to collect qc metrics. 
echo ====================
CMD="java -Xmx16g -jar $PICARD CollectMultipleMetrics \
                                    --INPUT $IN_BAM \
                                    --OUTPUT ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics \
                                    --REFERENCE_SEQUENCE $REFERENCE_GENOME \
                                    --PROGRAM null \
                                    --PROGRAM CollectSequencingArtifactMetrics \
                                    --PROGRAM CollectAlignmentSummaryMetrics \
                                    --PROGRAM CollectGcBiasMetrics \
                                    --PROGRAM QualityScoreDistribution \
                                    --VALIDATION_STRINGENCY LENIENT"
#CollectSequencingArtifactMetrics substitution_rate of the 6 base changes ,pre_adapter_summary_metrics TOTAL_QSCORE <35
#CollectAlignmentSummaryMetrics %PF_HQ_ALIGNED_READS(mapping quality of >20) %PF_HQ_ALIGNED_Q20_BASES PCT_READS_ALIGNED_IN_PAIRS PCT_PF_READS_IMPROPER_PAIRS PF_MISMATCH_RATE STRAND_BALANCE
#PF_NOISE_READS: number of PF reads that are composed entirly out of As and/or Ns
#PF_HQ_MEDIAN_MISMATCHES: median number of mismatches in PF_HQ_ALIGNED_READS
#PF_HQ_ERROR_RATE: percentage of bases that mismatch in PF_HQ_ALIGNED_READS
#BAD_CYCLES: Number of instrument cycles where 80% or more base calls were no-calls
#PCT_CHIMERAS: percentage of reads that map outside of a maximum insert size (usually 100kb) or that have the two ends mapping to different chromosomes
#PCT_ADAPTER: percentage of unaligned PF reads that match to a known adapter sequence right from the start of the read.
#CollectGcBiasMetrics AT_DROPOUT GC_DROPOUT

if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/multiple_metrics.alignment_summary_metrics || $REPLACE_TABLE != "false" ]];then
	echo $CMD
	eval $CMD
fi

CMD="java -Xmx16g -jar ${PICARD} QualityScoreDistribution \
                                  --INPUT ${IN_BAM} \
                                  --OUTPUT ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt \
                                  --CHART_OUTPUT ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.pdf \
                                  --VALIDATION_STRINGENCY LENIENT"
#QualityScoreDistribution can get the percentage of bases higher than 30
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt || $REPLACE_TABLE != "false" ]];then
	echo $CMD
	eval $CMD
fi

#Picard markduplicates has to be run before CollectHsMetrics to enable the HS_LIBRARY_SIZE field
CMD="java -Xmx16g -jar ${PICARD} MarkDuplicates \
                                  --INPUT ${IN_BAM} \
                                  --OUTPUT ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.bam \
                                  --METRICS_FILE ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.txt"
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.txt || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.txt"
  echo "CMD ->"
	echo $CMD
	eval $CMD
fi

#Collects hybrid-selection (HS) metrics for a SAM or BAM file. This tool takes a SAM/BAM file input and collects metrics that are specific for sequence datasets generated through hybrid-selection. Hybrid-selection (HS) is the most commonly used technique to capture exon-specific sequences for targeted sequencing experiments such as exome sequencing; for more information, please see the corresponding GATK Dictionary entry. 
#not very useful, similar function with coverage report
#AT_DROPOUT GC_DROPOUT
#FOLD_80_BASE_PENALTY: measure of non-uniformity. Number of units of sequencing necessary to raise 80% of the bases to the mean coverage (<5, between 3 and 4 for a good exome set)
#HS_LIBRARY_SIZE empty...
CMD="java -Xmx16g -jar ${PICARD} CollectWgsMetrics \
                                  --INPUT ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.bam \
                                  --OUTPUT ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt \
                                  --REFERENCE_SEQUENCE $REFERENCE_GENOME \
                                  --VALIDATION_STRINGENCY LENIENT"
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt"
  echo "CMD ->"
	echo $CMD
	eval $CMD
fi

#CollectOxoGMetrics #average OXIDATION_Q

CMD="java -Xmx16g -jar ${PICARD} CollectOxoGMetrics \
                                  --INPUT $IN_BAM \
                                  --OUTPUT ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt \
                                  --REFERENCE_SEQUENCE $REFERENCE_GENOME \
                                  --VALIDATION_STRINGENCY LENIENT"
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt"
  echo "CMD ->"
	echo $CMD
	eval $CMD
fi

CMD="java -Xmx16g -jar ${GATK} \
                        --analysis_type CallableLoci \
                        --reference_sequence $REFERENCE_GENOME \
                        --input_file $IN_BAM \
                        -summary ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt \
                        -o ${BUFFER_DIR}/PRE_QC/${NAME}/callable_status.bed"
#table.txt POOR_MAPPING_QUALITY (<10),LOW_COVERAGE(<4),NO_COVERAGE,still the callable_status.bed file can serve as a final vcf filter
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt"
  echo "CMD ->"
	echo $CMD
	eval $CMD
fi

CMD="java -Xmx16g -jar $GATK \
                        --analysis_type ErrorRatePerCycle \
                        --reference_sequence $REFERENCE_GENOME \
                        --input_file $IN_BAM \
                        -o ${BUFFER_DIR}/PRE_QC/${NAME}/error_rates.gatkreport.txt"
#not very useful, need to run across more samples to see the trend
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/error_rates.gatkreport.txt || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/error_rates.gatkreport.txt"
  echo "CMD ->"
	echo $CMD
	eval $CMD
fi

##filter in bam proper aligned reads first
#samtools view -@ 8 -f 3 -b $IN_BAM > ${TMP_DIR}/${NAME}_filtered.bam
#samtools index -@ 8 ${TMP_DIR}/${NAME}_filtered.bam

CMD="java -Xmx16g -jar $GATK \
                        --analysis_type CountTerminusEvent \
                        --reference_sequence $REFERENCE_GENOME \
                        --input_file ${TMP_DIR}/${NAME}_filtered.bam \
                        -o ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt "
#can try to add the two stats to the report
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt"
  echo "CMD ->"

  #filter in bam proper aligned reads first (move samtools command line from outside to inside "if")
  samtools view -@ 8 -f 3 -b $IN_BAM > ${TMP_DIR}/${NAME}_filtered.bam
  samtools index -@ 8 ${TMP_DIR}/${NAME}_filtered.bam

	echo $CMD
	eval $CMD

	rm -rf ${TMP_DIR}/${NAME}_filtered.bam* 2> /dev/null
fi
#rm -rf ${TMP_DIR}/${NAME}_filtered.bam* 2> /dev/null

CMD="java -Xmx16g -jar $GATK \
                        -T ReadLengthDistribution \
                        -R $REFERENCE_GENOME \
                        -I $IN_BAM \
                        -o ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl"
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl"
  echo "CMD ->"
	echo $CMD
	eval $CMD
fi

CMD="java -Xmx16g -jar $PICARD CollectInsertSizeMetrics \
                                --INPUT $IN_BAM \
                                --OUTPUT ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt \
                                --Histogram_FILE ${BUFFER_DIR}/PRE_QC/${NAME}/InsertSizeHist.pdf \
                                --DEVIATIONS 10.0 \
                                --MINIMUM_PCT 0.05 \
                                --VALIDATION_STRINGENCY LENIENT"
if [[ ! -f ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt || $REPLACE_TABLE != "false" ]];then
  echo
  echo "Output: ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt"
  echo "CMD ->"
	echo $CMD
	eval $CMD
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
TOTAL_BASES=`awk 'NR>8&&NR<16' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt | awk -F "\t" '{sum+=$2} END {print sum}'`
#BASES_Q_AVE=`awk -F"\t" '{if(NR>8 && NR<59){sum1+=$2;sum2+=$1*$2}} END{print sum2/sum1}' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt | awk -F "\t" '{sum+=$2} END {print sum}'`
BASES_Q_AVE=`awk -F '\t' '{if(NR>8 && NR<59){sum1+=$2;sum2+=$1*$2}} END{print sum2/sum1}' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt`
#BASES_Q_12DOWN=`awk 'NR>8&&NR<11' ${BUFFER_DIR}/PRE_QC/${NAME}/qual_score_dist.txt | awk -F "\t" '{sum+=$2} END {print sum}'`
if [[ ! $TOTAL_BASES -gt 0 ]]; then	
	echo "Error. TOTAL_BASES is less than 0, something is wrong."
#else
#	PERCENT_BASES_Q_32UP=`echo "$BASES_Q_32UP / $TOTAL_BASES" | bc -l | xargs -I {} printf "%5.4f" {}`
#	PERCENT_BASES_Q_12DOWN=`echo "$BASES_Q_12DOWN / $TOTAL_BASES" | bc -l | xargs -I {} printf "%5.4f" {}`
fi

#Collects hybrid-selection (HS) metrics AT_DROPOUT GC_DROPOUT FOLD_80_BASE_PENALTY HS_LIBRARY_SIZE
MEAN_COVERAGE=`awk -F "\t" '{if(NR==8){print $2}}' ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt`
SD_COVERAGE=`awk -F "\t" '{if(NR==8){print $3}}' ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt`
MEDIAN_COVERAGE=`awk -F "\t" '{if(NR==8){print $4}}' ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt`
PCT_EXC_MAPQ=`awk -F "\t" '{if(NR==8){print $6}}' ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt`
PCT_EXC_DUPE=`awk -F "\t" '{if(NR==8){print $7}}' ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt`
PCT_EXC_UNPAIRED=`awk -F "\t" '{if(NR==8){print $8}}' ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt`
PCT_EXC_OVERLAP=`awk -F "\t" '{if(NR==8){print $10}}' ${BUFFER_DIR}/PRE_QC/${NAME}/wgs_metrics.txt`

#Collects Estimated_library_size from markduplicates
ESTIMATED_LIBRARY_SIZE=`head -n 8 ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.txt | tail -n 1| awk -F "\t" '{print $10}'`

#CollectOxoGMetrics average OXIDATION_Q

OXIDATION_Q_AVE=`awk 'NR>7&&NR<24' ${BUFFER_DIR}/PRE_QC/${NAME}/oxoG_metrics.txt |  awk '{ sum += $12 } END { print sum /16 }'`

#GATK CallableLoci  POOR_MAPPING_QUALITY (<10),LOW_COVERAGE(<4),NO_COVERAGE
LOW_COVERAGE_CAPTURE_REGION=`awk '{if (NR==5){print $2}}' ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt`
NO_COVERAGE_CAPTURE_REGION=`awk '{if (NR==4){print $2}}' ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt`
POOR_MAPPING_QUALITY=`awk '{if (NR==7){print $2}}' ${BUFFER_DIR}/PRE_QC/${NAME}/table.txt`

PERCENT_LOW_COVERAGE_CAPTURE_REGION=`echo "$LOW_COVERAGE_CAPTURE_REGION / $TOTAL_BASE " | bc -l | xargs -I {} printf "%5.4f" {}`
PERCENT_NO_COVERAGE_CAPTURE_REGION=`echo "$NO_COVERAGE_CAPTURE_REGION / $TOTAL_BASE " | bc -l | xargs -I {} printf "%5.4f" {}`
PERCENT_POOR_MAPPING_QUALITY=`echo "$POOR_MAPPING_QUALITY / $TOTAL_BASE " | bc -l | xargs -I {} printf "%5.4f" {}`

#GATK CountTerminusEvent 
READS_ENDING_IN_INDELS=`head -1 ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt | awk -F ": " '{print $2}'`
READS_ENDING_IN_SOFTCLIPS=`head -2 ${BUFFER_DIR}/PRE_QC/${NAME}/output.txt  | tail -1 | awk -F ": " '{print $2}'`

#GATK ReadLengthDistribution
TOTAL_READS=`awk 'NR>4' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl | awk '{sum+=$NF} END {print sum}'`
TOTAL_BASES=`awk 'NR>4' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl | awk '{sum+=$2*$1} END {print sum}'`
TOTAL_READS_WITH_3LESS_CUT=`awk '/^$/ {nlstack=nlstack "\n";next;} {printf "%s",nlstack; nlstack=""; print;}' ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_rld.tbl | tail -4 | awk '{sum+=$2} END {print sum}'`
if [[ ! $TOTAL_READS -gt 0 ]]; then	
	echo "Error. TOTAL_READS is less than 0, something is wrong."
else
	AVERAGE_READ_LENGTH=`echo "$TOTAL_BASES / $TOTAL_READS" | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_TOTAL_READS_WITH_3LESS_CUT=`echo "$TOTAL_READS_WITH_3LESS_CUT / $TOTAL_READS " | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_READS_ENDING_IN_INDELS=`echo "$READS_ENDING_IN_INDELS / $TOTAL_READS " | bc -l | xargs -I {} printf "%5.4f" {}`
	PERCENT_READS_ENDING_IN_SOFTCLIPS=`echo "$READS_ENDING_IN_SOFTCLIPS / $TOTAL_READS " | bc -l | xargs -I {} printf "%5.4f" {}`
fi

MEDIAN_INSERT_SIZES=`awk -F '\t' '{ if (NR==8) print $1 }' ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt`
MEAN_INSERT_SIZES=`awk -F '\t' '{ if (NR==8) print $6 }' ${BUFFER_DIR}/PRE_QC/${NAME}/insertsize_metrics.txt`

L="${NAME}\t${AVERAGE_READ_LENGTH}\t${PERCENT_TOTAL_READS_WITH_3LESS_CUT}\t${PF_HQ_MEDIAN_MISMATCHES}\t${PF_HQ_MISMATCH_RATE_F}\t${PF_HQ_MISMATCH_RATE_R}\t${STRAND_BALANCE_F}\t${STRAND_BALANCE_R}\t${BAD_CYCLES}\t${PCT_CHIMERAS}\t${PCT_ADAPTER}"\
"\t${PERCENT_PF_HQ_ALIGNED_READS}\t${PERCENT_PF_HQ_ALIGNED_BASES_L}\t${PERCENT_PF_HQ_ALIGNED_Q20_BASES_L}\t${PERCENT_PF_HQ_ALIGNED_BASES_R}\t${PERCENT_PF_HQ_ALIGNED_Q20_BASES_R}\t${BASES_Q_AVE}"\
"\t${AT_DROPOUT}\t${GC_DROPOUT}\t${MEAN_COVERAGE}\t${SD_COVERAGE}\t${MEDIAN_COVERAGE}\t${PCT_EXC_MAPQ}\t${PCT_EXC_UNPAIRED}\t${PCT_EXC_DUPE}\t${PCT_EXC_OVERLAP}\t${ESTIMATED_LIBRARY_SIZE}"\
"\t${LOWEST_PREADATPER_TOTAL_QSCORE_BASE}\t${LOWEST_PREADATPER_TOTAL_QSCORE}\t${LOWEST_BAITBIAS_TOTAL_QSCORE_BASE}\t${LOWEST_BAITBIAS_TOTAL_QSCORE}\t${OXIDATION_Q_AVE}"\
"\t${PERCENT_LOW_COVERAGE_CAPTURE_REGION}\t${PERCENT_NO_COVERAGE_CAPTURE_REGION}\t${PERCENT_POOR_MAPPING_QUALITY}\t${PERCENT_READS_ENDING_IN_INDELS}\t${PERCENT_READS_ENDING_IN_SOFTCLIPS}\t${MEDIAN_INSERT_SIZES}\t${MEAN_INSERT_SIZES}"

echo "Contents to be streamed to the flowcell-level report file:"
echo -e $L
#echo -e $L >> $OUT_REPORT

#remove all intermediate bam file to save space
strFile="${TMP_DIR}/${NAME}_filtered.bam"
if [ -f ${strFile} ]; then
  rm ${strFile}
fi

strFile="${TMP_DIR}/${NAME}_filtered.bam.bai"
if [ -f ${strFile} ]; then
  rm ${strFile}
fi

strFile="${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.bam"
if [ -f ${strFile} ]; then
  rm ${strFile}
fi

strFile="${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.bai"
if [ -f ${strFile} ]; then
  rm ${strFile}
fi

#rm ${TMP_DIR}/${NAME}_filtered.bam ${TMP_DIR}/${NAME}_filtered.bam.bai ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.bam ${BUFFER_DIR}/PRE_QC/${NAME}/${NAME}_dedup.bai
echo ====================
date
echo Done!
echo ====================

# update flag
rm ${strFlagWorking}
touch ${strFlagDone}
