#!/bin/bash
set -o pipefail
# Note that this version turn off set flags for testing purpose
#
# 
# map pair-ended, Illumina 1.8+ format fastq files to reference genome, generate *_original.bam
# then remove optical and PCR duplicates, generate *_dedup.bam
# and further remove incorrectly mapped pairs, generate *_dedup_correctly_paired.bam
# then lastly generate the flagstat report of these three bam files for later report generation purpose
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
# calls NovoAlign to align paired-end reads to reference genome
# PATH=$NOVOCRAFT_PRG_DIR:$PATH

#takes in arguments
argIndex=1
for arg in "$@"
do
  if [[ ${argIndex} == 1 ]]; then
    FASTQ1=${arg}
  elif [[ ${argIndex} == 2 ]]; then
    FASTQ2=${arg}
  elif [[ ${argIndex} == 3 ]]; then
    OUTBAM=${arg}
  elif [[ ${argIndex} == 4 ]]; then
    RUN_ID=${arg}
  elif [[ ${argIndex} == 5 ]]; then
    sampleType=${arg}
  elif [[ ${argIndex} == 6 ]]; then
    TmpFolder=${arg}
  elif [[ ${argIndex} == 7 ]]; then
    JOB_NAME=${arg}
  elif [[ ${argIndex} == 8 ]]; then
    reference=${arg}
  elif [[ ${argIndex} == 9 ]]; then
    refSamIndex=${arg}
  elif [[ ${argIndex} == 10 ]]; then
    refSamBed=${arg}
  elif [[ ${argIndex} == 11 ]]; then
    runFlowcellDir=${arg}
  elif [[ ${argIndex} == 12 ]]; then
    alignerTool=${arg}
  elif [[ ${argIndex} == 13 ]]; then
    refVersion=${arg}
  elif [[ ${argIndex} == 14 ]]; then
    sizeRam=${arg}
  elif [[ ${argIndex} == 15 ]]; then
    doneAlignment=${arg}
  elif [[ ${argIndex} == 16 ]]; then
    doneDedup=${arg}
  elif [[ ${argIndex} == 17 ]]; then
    donePhix=${arg}
  elif [[ ${argIndex} == 18 ]]; then
    doneUMI=${arg}
  elif [[ ${argIndex} == 19 ]]; then
    allLaneDemultiplexOnly=${arg}
  elif [[ ${argIndex} == 20 ]]; then
    coreNum=${arg}
  elif [[ ${argIndex} == 21 ]]; then
    flagWorking=${arg}
  elif [[ ${argIndex} == 22 ]]; then
    flagDone=${arg}
  elif [[ ${argIndex} == 22 ]]; then
    iMergedSample=${arg}  
  fi
  argIndex=$((argIndex + 1))
done

ILLUMINA_FLAG_SUB_DIR="Flag"
# The 3rd argument is supposed to pass only the prefix of the bam file including directory name
# e.g. "/CGF/Sequencing/Illumina/GA_PostRun_Analysis/GA_Data/121120_SN7001343_0058_AC1DV8ACXX/BAM/SB808862_GTGAAA_L008"
# but just in case the full name including extension is passed, strip down the extension
OUTDIR=`dirname $OUTBAM`
OUTNAME=`basename $OUTBAM .bam`
echo "OUTDIR : "${OUTDIR}
echo "OUTNAME: "${OUTNAME}
# the first token delimited by "_" is assumed to be the name of the sample
SAMPLE=`echo $OUTNAME |rev|cut -d_ -f5-|rev`
#remove the last 3part lane (https://unix.stackexchange.com/questions/234432/how-to-delete-the-last-column-of-a-file-in-linux)
LIB=`echo $OUTNAME | rev|cut -d_ -f4-|rev`

FLOWCELL_Folder_Name=`basename ${runFlowcellDir}`
FLOWCELL=`echo ${FLOWCELL_Folder_Name} | awk -F '_' '{print $NF}'`
FLOWCELL_DATE=`echo ${FLOWCELL_Folder_Name} | awk -F '_' '{print $1}'`

LANE=`echo $OUTNAME | sed 's/_HQ_paired//' |  tail -c 2`
INDEX=`echo $LIB | awk -F "_" '{print $NF}'`
# these are the flag files to indicate the process status
INDIR=`dirname $FASTQ1`


#RUN_ID=`echo $INDIR | cut -d/ -f8`
flagDir=${INDIR}

strSuffixDir=""

if [[ ${allLaneDemultiplexOnly} == "yes" ]]; then
  strSuffixDir="${DEMULTIPLEX_ONLY_SUB_DIR}"
elif [[ ${alignerTool} != "" ]] && [[ ${refVersion} != "" ]]; then
  strSuffixDir="${alignerTool}/${refVersion}"
fi

if [[ ${strSuffixDir} != "" ]]; then
  flagDir=${flagDir}/${ILLUMINA_FLAG_SUB_DIR}/${strSuffixDir}
fi

#if ! [[ -d ${flagDir} ]]; then
#  mkdir -p ${flagDir}
#fi

TMP_DIR=${ILLUMINA_PROCESSED_DATA_ROOT_DIR}/${RUN_ID}/ttemp
if [[ ${strSuffixDir} != "" ]]; then
  TMP_DIR=${TMP_DIR}/${strSuffixDir}
fi
if ! [[ -d ${TMP_DIR} ]]; then
  mkdir -p ${TMP_DIR}
fi

#Check if Current Sample Folder contains UMI Fastq File -->
bContainUMIReads=false
fastqFile=${INDIR}/*_UMI.fastq.gz
umiFastq=$(find ${INDIR} -type f -name "*_UMI.fastq.gz" | head -1)
if [[ -f ${umiFastq} ]]; then
  bContainUMIReads=true
fi
echo "INDIR   : "${INDIR}
echo "umiFastq: "${umiFastq}
#<--

# call novoalign program to generate *_original.bam
# per lab, the average insertion size of the PE is about 300bp
# allowing 100bp window -- "-i PE 300,100"
OUTBAM_ORIGINAL=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_ORIGINAL}.bam
OUTBAI_ORIGINAL=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_ORIGINAL}.bai

#determine the reference file based on sample is whole exome sample or whole genome sample
if [[ ${sampleType} != "DNAWholeGenome" ]] && [[ ${sampleType} != "DNAExome" ]]; then
  echo "Error: unkown sample type: ${sampleType}!"
  exit 1
fi

if [[ ! -f $OUTBAM_ORIGINAL ]] || [[ ! -f $OUTBAI_ORIGINAL ]]; then
  echo "===================="
  date
  echo "$(date) Reference genome mapping based on NovoAlign"
  echo "Input FASTQs: $FASTQ1 $FASTQ2"
  echo "Output BAM: $OUTBAM_ORIGINAL"
  echo "===================="

  # Step 1: ----- Run Novoalign ---->
  #1: Define parameters
  readsType="PE"
  format="ILM1.8"
  fragmentLen=300
  stdDeviation=125
  outputType="SAM"
  RGArea="@RG\tCN:CGR\tPL:ILLUMINA\tID:${FLOWCELL}.${LANE}\tSM:${SAMPLE}\tPU:${FLOWCELL}20${FLOWCELL_DATE}.${LANE}.${INDEX}"
  samFile="${TMP_DIR}/${OUTNAME}.sam"
  tmpRawBAM="${TMP_DIR}/${OUTNAME}.bam"
  #2: print version of novoalign
  #3: run novoalign
  if [[ ${doneAlignment} == "no" ]] || [[ ! -f ${samFile} ]]; then
    #Check Different Aligner
    if [[ ${alignerTool} == "NovoAlign" ]]; then
      echo "Run novoalign -->"
      echo "NovoAlign Version: "`${NOVOALIGN} --version`
      CMD="${NOVOALIGN} -c ${coreNum} \
                        -d ${reference} \
                        -f ${FASTQ1} ${FASTQ2} \
                        -F ${format} \
                        -i ${readsType} ${fragmentLen},${stdDeviation} \
                        -o ${outputType} $'${RGArea}' > ${samFile}"
      echo ${CMD}
      SECONDS=0
      eval ${CMD}
      result=$?
      duration=$SECONDS
      echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
      if [[ ${result} -ne 0 ]]; then
        echo "Error: $(date) NovoAlign was failed!"
        exit 1
      else
        echo "$(date) NovoAlign alignment was finished successfully!"
      fi
    elif [[ ${alignerTool} == "BWA" ]]; then
      echo "Do Alignment by using BWA"
      CMD="bwa mem -t ${coreNum} -R \"${RGArea}\" ${reference} ${FASTQ1} ${FASTQ2} > ${samFile}"
      echo ${CMD}
      SECONDS=0
      eval ${CMD}
      result=$?
      duration=$SECONDS
      echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
      if [[ ${result} -ne 0 ]]; then
        echo "Error: $(date) BWA was failed!"
        exit 1
      else
        echo "$(date) BWA alignment was finished successfully!"
      fi
    else
      echo "Error: Aligner is invalid!"
      exit 1
    fi
  else
    echo "$(date) ${alignerTool} alignment was finished successfully!"
  fi
  #<----------------

  #Step 2: sam converting to bam
  tmpSortedBAM="${TMP_DIR}/${OUTNAME}_sorted.bam"
  echo
  echo "$(date) sam converting to bam"
  if [[ ! -f ${tmpSortedBAM} ]]; then
    CMD="samtools view -bt ${refSamIndex} -o ${tmpRawBAM} ${samFile}"
    echo $CMD
    SECONDS=0
    eval $CMD
    result=$?
    duration=$SECONDS
    echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
    if [[ ${result} -ne 0 ]]; then
      echo "Error: $(date) sam converting to bam was failed!"
      exit 1
    else
      echo "$(date) sam converting to bam was finished successfully!"
    fi
  else
    echo "$(date) sam converting to bam was finished successfully!"
  fi

  #Step 3: bam to sorted bam
  tmpRenamedBAM="${TMP_DIR}/${OUTNAME}_renamed.bam"
  echo
  echo "$(date) samtools sorting bam"
  if [[ ! -f ${tmpRenamedBAM} ]]; then
    CMD="samtools sort -@ ${coreNum} -T ${TMP_DIR}/${OUTNAME}_temp -o ${tmpSortedBAM} ${tmpRawBAM}"
    echo $CMD
    SECONDS=0
    eval $CMD
    result=$?
    duration=$SECONDS
    echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
    if [[ ${result} -eq 0 ]]; then
      echo "$(date) samtools sorting bam was finished successfully!"
    else
      echo ${result}
      echo "Error: $(date) samtools sorting bam was failed!"
      exit 1;
    fi
  else
    echo "$(date) samtools sorting bam was finished successfully!"
  fi

  #Step 4: Refine RG Group in Sorted Bam
  echo
  CMD="$PICARD_ADDORREPLACEREADGROUPS \
        --INPUT ${TMP_DIR}/${OUTNAME}_sorted.bam \
        --OUTPUT ${tmpRenamedBAM} \
        --RGID "${FLOWCELL}.${LANE}" \
        --RGLB "${LIB}" \
        --RGSM "${SAMPLE}" \
        --RGPL "ILLUMINA" \
        --RGPU "${FLOWCELL}20${FLOWCELL_DATE}.${LANE}.${INDEX}" \
        --RGCN "CGR" \
        --SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY LENIENT \
        --CREATE_INDEX true"
  echo $CMD
  SECONDS=0
  eval $CMD
  result=$?
  duration=$SECONDS
  echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
  if [[ ${result} -ne 0 ]]; then
    echo "Error: $(date) Picard AddOrReplaceReadGroups was failed!"
    exit 1
  else
      echo "$(date) Picard AddOrReplaceReadGroups was finished successfully!"
  fi

  #Step 5: change file name
  echo
  CMD="mv ${TMP_DIR}/${OUTNAME}_renamed.bam ${OUTBAM_ORIGINAL} &&
       mv ${TMP_DIR}/${OUTNAME}_renamed.bai ${OUTBAI_ORIGINAL}"
  echo $CMD
  eval $CMD
  if [[ $? -ne 0 ]]; then
    echo "Error: $(date) sam indexing was failed!"
    exit 1
  else
      echo "$(date) sam indexing was finished successfully!"
  fi
  CMD="rm ${TMP_DIR}/${OUTNAME}*"
  echo $CMD
  eval $CMD

  if [[ ! -f $OUTBAM_ORIGINAL ]] ; then
    echo "$OUTBAM_ORIGINAL file expected as result but not found. Novoalign runtime exception?" 1>&2
    exit 1
  fi
else
#  echo "Error: $OUTBAM_ORIGINAL or $OUTBAI_ORIGINAL is alreay existed!"
#  exit 1;
  echo "Do Next Step: both $OUTBAM_ORIGINAL and $OUTBAI_ORIGINAL are alreay existed!"
fi

# call picard to remove optical and PCR duplicates
OUTBAM_DEDUP=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_DEDUP}.bam
OUTBAI_DEDUP=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_DEDUP}.bai
echo ====================
date
echo $(date) Removing optical and PCR duplicates
echo Input BAM: $OUTBAM_ORIGINAL
echo Output BAM: $OUTBAM_DEDUP
echo ====================

#Get Metrics File Name --->
reportDir=${runFlowcellDir}/${ILLUMINA_REPORT_SUB_DIR}
if [[ ${strSuffixDir} != "" ]]; then
  reportDir=${reportDir}/${strSuffixDir}
fi
if ! [[ -d ${reportDir} ]]; then
  mkdir -p ${reportDir}
fi
dupMetricsFile=${reportDir}/${OUTNAME}_duplication_metrics.txt
#<---
FLOWCELL_SERIAL=`echo $RUN_ID|cut -d'_' -f2`

if [[ ${doneDedup} == "no" ]] || [[ ! -f ${OUTBAM_DEDUP} ]]; then
  SECONDS=0
  if [[ $FLOWCELL_SERIAL = "K00278" ]] || [[ $FLOWCELL_SERIAL = "A00423" ]]; then
    echo "It's a patterned flow cell"
    CMD="$PICARD_MARK_DUPLICATES \
          --INPUT $OUTBAM_ORIGINAL \
          --OUTPUT $OUTBAM_DEDUP \
          --METRICS_FILE ${dupMetricsFile} \
          --REMOVE_DUPLICATES true \
          --CREATE_INDEX true \
          --ASSUME_SORTED true \
          --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500  \
          --VALIDATION_STRINGENCY STRICT \
          --TMP_DIR ${TmpFolder}"
    echo $CMD
    eval $CMD
  else
    echo "It's a non patterned flow cell"
    CMD="$PICARD_MARK_DUPLICATES \
          --INPUT $OUTBAM_ORIGINAL \
          --OUTPUT $OUTBAM_DEDUP \
          --METRICS_FILE ${dupMetricsFile} \
          --REMOVE_DUPLICATES true \
          --CREATE_INDEX true \
          --ASSUME_SORTED true \
          --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
          --VALIDATION_STRINGENCY STRICT \
          --TMP_DIR $TmpFolder"
    echo $CMD
    eval $CMD
  fi
  result=$?
  duration=$SECONDS
  echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

  if [[ ${result} -eq 0 ]]; then
    echo "$(date) Dedup was finished successfully!"
  else
    echo "Error: $(date) dedup was failed!"
    exit 1
  fi

  if [[ ! -f ${OUTBAM_DEDUP} ]]; then
    echo "$OUTBAM_DEDUP file expected as result but not found. Picard runtime exception?" 1>&2
    exit 1
  fi
else
  if [[ ! -f ${OUTBAM_DEDUP} ]]; then
    echo "$OUTBAM_DEDUP file expected as result but not found. Picard runtime exception?" 1>&2
    exit 1
  else
    echo "$(date) Dedup was finished successfully!"
  fi
fi

echo ====================
date
echo Calculating Insertzises for dedup BAM
echo Input BAM: $OUTBAM_DEDUP
echo Output FILE: ${OUTNAME}_insertsize_metrics.txt
echo ====================
#Get Metrics File Name --->
insertMetricsFile=${reportDir}/${OUTNAME}_insertsize_metrics.txt
histogramFile=${reportDir}/${OUTNAME}_InsertSizeHist.pdf
#<---
CMD="$PICARD_INSERTSIZE \
      --INPUT $OUTBAM_DEDUP \
      --OUTPUT ${insertMetricsFile} \
      --Histogram_FILE ${histogramFile} \
      --DEVIATIONS 10.0 \
      --MINIMUM_PCT 0.05 \
      --VALIDATION_STRINGENCY STRICT"
echo $CMD
SECONDS=0
eval $CMD
result=$?
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

if [[ ${result} -eq 0 ]]; then
  echo "$(date): Insertsizes calculation for dedup BAM was finished successfully!"
else
  echo "Error: $(date) insertsizes for dedup BAM was failed!"
  exit 1
fi

# remove phiX control records
OUTBAM_NOPHIX=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_DEDUP}_${SUFFIX_BAM_NOPHIX}.bam
OUTBAI_NOPHIX=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_DEDUP}_${SUFFIX_BAM_NOPHIX}.bai
echo ====================
date
echo Remove phiX records
echo Input BAM: $OUTBAM_DEDUP
echo Output BAM: $OUTBAM_NOPHIX
echo ====================

#if [[ $hasWGS -gt 0 ]]; then
if [[ ${sampleType} == "DNAWholeGenome" ]]; then
  echo "Whole genome data, the reference file for whole genome data does not have phiX contig, no need to remove phiX reads, skip removing phiX records step and rename BAM file."
  mv $OUTBAM_DEDUP $OUTBAM_NOPHIX
  mv $OUTBAI_DEDUP $OUTBAI_NOPHIX
  if [[ $? -eq 0 ]]; then
    echo "Rename $OUTBAM_DEDUP to $OUTBAM_NOPHIX is successful!"
	echo "Rename $OUTBAI_DEDUP to $OUTBAI_NOPHIX is successful!"
  else
    echo "Error: Rename $OUTBAM_DEDUP to $OUTBAM_NOPHIX or rename $OUTBAI_DEDUP to $OUTBAI_NOPHIX is failed!"
    exit 1;
  fi
else
  if [[ ${donePhix} == "no" ]] || [[ ! -f ${OUTBAM_NOPHIX} ]]; then
    CMD="samtools view -h -L ${refSamBed} -U ${OUTDIR}/${OUTNAME}_phiX174.bam $OUTBAM_DEDUP \
          | awk '{if ((\$2 !~ /$PHIX_REFERENCE_ID_TO_CLEAN/) && (\$7 !~ /$PHIX_REFERENCE_ID_TO_CLEAN/))  print}' \
          | samtools view -bS -o $OUTBAM_NOPHIX - && samtools index $OUTBAM_NOPHIX $OUTBAI_NOPHIX"
    echo $CMD
    SECONDS=0
    eval $CMD
    result=$?
    duration=$SECONDS
    echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
    if [[ ${result} -eq 0 ]]; then
      echo "Phix records were removed successfully!"
    else
      echo "Error: Phix records removal was failed!"
      exit 1;
    fi
    rm ${OUTDIR}/${OUTNAME}_phiX174.bam
    if [[ ! -f $OUTBAM_NOPHIX ]]; then
      echo "$OUTBAM_NOPHIX file expected as result but not found. Samtools runtime exception?" 1>&2
      exit 1
    fi
  else
    if [[ ! -f $OUTBAM_NOPHIX ]]; then
      echo "$OUTBAM_NOPHIX file expected as result but not found. Samtools runtime exception?" 1>&2
      exit 1
    else
      echo "Phix records were removed successfully!"
    fi
  fi
fi

# generate flagstat reports for different levels of BAM files
FLAGSTAT_FILE=${reportDir}/${OUTNAME}_flagstats.txt
echo ====================
date
echo Generating BAM flagstat reports
echo Output: $FLAGSTAT_FILE
echo ====================

echo $OUTBAM_ORIGINAL > $FLAGSTAT_FILE
samtools flagstat -@ ${coreNum} $OUTBAM_ORIGINAL >> $FLAGSTAT_FILE
result=$?
if [[ ${result} -ne 0 ]]; then
  echo "Error: samtools flagstat(original) was failed!"
  exit 1;
fi
echo " " >> $FLAGSTAT_FILE
echo " " >> $FLAGSTAT_FILE

#if [[ $hasWGS -gt 0 ]]; then
if [[ ${sampleType} == "DNAWholeGenome" ]]; then
  echo $OUTBAM_NOPHIX >> $FLAGSTAT_FILE
  samtools flagstat -@ ${coreNum} $OUTBAM_NOPHIX >> $FLAGSTAT_FILE
  result=$?
  if [[ ${result} -ne 0 ]]; then
    echo "Error: samtools flagstat(NOPHIX) was failed!"
    exit 1;
  fi
  echo " " >> $FLAGSTAT_FILE
  echo " " >> $FLAGSTAT_FILE
else
  echo $OUTBAM_DEDUP >> $FLAGSTAT_FILE
  samtools flagstat -@ ${coreNum} $OUTBAM_DEDUP >> $FLAGSTAT_FILE
  result=$?
  if [[ ${result} -ne 0 ]]; then
    echo "Error: samtools flagstat(DEDUP) was failed!"
    exit 1;
  fi
  echo " " >> $FLAGSTAT_FILE
  echo " " >> $FLAGSTAT_FILE
fi

echo $OUTBAM_NOPHIX >> $FLAGSTAT_FILE
samtools flagstat -@ ${coreNum} $OUTBAM_NOPHIX >> $FLAGSTAT_FILE
result=$?
if [[ ${result} -ne 0 ]]; then
  echo "Error: samtools flagstat(NOPHIX_second) was failed!"
  exit 1;
fi
echo " " >> $FLAGSTAT_FILE
echo " " >> $FLAGSTAT_FILE

if [[ ! -f $FLAGSTAT_FILE ]] ; then
  echo $FLAGSTAT_FILE file expected as result but not found. Samtools runtime exception? 1>&2
  exit 1
fi

#-> Calculate number of mapped reads for both HQ Mapped Bam and Dedupe Bam respectively (Not number of Alignment)
#(1) For HQ Mapped BAM
echo
echo "===================="
date
echo "Count Mapped Reads for HQ Mapped BAM"
echo "Output: $FLAGSTAT_FILE"
echo "===================="
SECONDS=0
#left mate
tmpBAM=${OUTBAM_ORIGINAL}
mappedNumLeftMate=$(samtools view -@ ${coreNum} -f 0x40 -F 0x4 ${tmpBAM} | cut -f1 | sort | uniq | wc -l)
#right mate
mappedNumRightMate=$(samtools view -@ ${coreNum} -f 0x80 -F 0x4 ${tmpBAM} | cut -f1 | sort | uniq  | wc -l)
mappedNumTotal=$((${mappedNumLeftMate} + ${mappedNumRightMate}))
echo "Total Number of Mapped Reads in Original BAM:${mappedNumTotal}" >> $FLAGSTAT_FILE
echo " " >> $FLAGSTAT_FILE
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

#(2) For HQ Dedup BAM
echo
echo "===================="
date
echo "Count Mapped Reads for Dedup BAM"
echo "Output: $FLAGSTAT_FILE"
echo "===================="
SECONDS=0
tmpBAM=${OUTBAM_DEDUP}
if [[ ${sampleType} == "DNAWholeGenome" ]]; then
  tmpBAM=${OUTBAM_NOPHIX}
fi
#left mate
mappedNumLeftMate=$(samtools view -@ ${coreNum} -f 0x40 -F 0x4 ${tmpBAM} | cut -f1 | sort | uniq | wc -l)
#right mate
mappedNumRightMate=$(samtools view -@ ${coreNum} -f 0x80 -F 0x4 ${tmpBAM} | cut -f1 | sort | uniq  | wc -l)
mappedNumTotal=$((${mappedNumLeftMate} + ${mappedNumRightMate}))
echo "Total Number of Mapped Reads in Dedup BAM:${mappedNumTotal}" >> $FLAGSTAT_FILE
echo " " >> $FLAGSTAT_FILE
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
#<--

#(3) For DedupNoPhix BAM (only do this for the case of NONE WGS flowcell)
if ! [[ ${sampleType} == "DNAWholeGenome" ]]; then
  echo
  echo "===================="
  date
  echo "Count Mapped Reads for DedupNoPhix BAM"
  echo "Output: $FLAGSTAT_FILE"
  echo "===================="
  SECONDS=0
  tmpBAM=${OUTBAM_NOPHIX}
  #left mate
  mappedNumLeftMate=$(samtools view -@ ${coreNum} -f 0x40 -F 0x4 ${tmpBAM} | cut -f1 | sort | uniq | wc -l)
  #right mate
  mappedNumRightMate=$(samtools view -@ ${coreNum} -f 0x80 -F 0x4 ${tmpBAM} | cut -f1 | sort | uniq  | wc -l)
  mappedNumTotal=$((${mappedNumLeftMate} + ${mappedNumRightMate}))
  echo "Total Number of Mapped Reads in DedupNoPhix BAM:${mappedNumTotal}" >> $FLAGSTAT_FILE
  echo " " >> $FLAGSTAT_FILE
  duration=$SECONDS
  echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
fi
#<---

# delete intermediate files to save space
if [ $NOVOALIGN_DELETE_INTERMEDIATE_FILES = true ]; then
  rm -rf $OUTBAM_DEDUP $OUTBAI_DEDUP
fi

#--> Add UMI Code and save it as "dedup" + "npphix" + "UMI" (the softlink also point to this file later)
if ${bContainUMIReads}; then
  OUTBAM_DEDUP_NOPHIX_UMI=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_DEDUP}_${SUFFIX_BAM_NOPHIX}_UMI.bam
  OUTBAI_DEDUP_NOPHIX_UMI=${OUTDIR}/${OUTNAME}_${SUFFIX_BAM_DEDUP}_${SUFFIX_BAM_NOPHIX}_UMI.bai
  echo ====================
  date
  echo Add UMI Flag
  echo Input BAM: ${OUTBAM_NOPHIX}
  echo Output BAM: ${OUTBAM_DEDUP_NOPHIX_UMI}
  echo ====================

  #Estimate the size of RAM
  FGBio="/DCEG/CGF/Bioinformatics/software/lix33/fgbio/fgbio-1.2.0.jar"
  if [[ ${doneUMI} == "no" ]] || [[ ! -f ${OUTBAM_DEDUP_NOPHIX_UMI} ]]; then
    #Run FgBio
    SECONDS=0
    echo "Do UMI Annotation (Notice we should actually NEED the option '-Xmx')->"
    CMD="java -Xmx${sizeRam} -jar ${FGBio} AnnotateBamWithUmis \
                             -i ${OUTBAM_NOPHIX} \
                             -f ${umiFastq} \
                             -o ${OUTBAM_DEDUP_NOPHIX_UMI}"
    echo ${CMD}
    eval ${CMD}
    result=$?
    #print running time
    duration=$SECONDS
    echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

    if [[ ${result} -eq 0 ]]; then
      echo "UMI Annotation has finished Successfully!"
    else
      echo "Error: UMI Annotation was Failed!"
      exit 1
    fi
  else
    if [[ ! -f ${OUTBAM_DEDUP_NOPHIX_UMI} ]]; then
      echo "${OUTBAM_DEDUP_NOPHIX_UMI} file expected as result but not found. FGBio runtime exception?" 1>&2
      exit 1
    else
      echo "UMI Annotation has finished Successfully!"
    fi
  fi

#  #Copy Bai --> We do not need to copy it, since "AnnotateBamWithUmis" will generate both
#  cp ${OUTBAI_NOPHIX} ${OUTBAI_DEDUP_NOPHIX_UMI}
#  #<--
  #Notice: we still keep the original OUTBAM_NOPHIX and its BAI file in this version
  echo
fi
#<--

# create soft links to the two file for compatibility to legacy scripts
OUTNAME_SIMPLE=`echo $OUTNAME | rev | cut -d_ -f3- | rev`

softLinkFolder=${OUTDIR}
softLinkBAM=${OUTNAME_SIMPLE}.bam
softLinkBAI=${OUTNAME_SIMPLE}.bai

rm -f ${softLinkFolder}/${softLinkBAM}
rm -f ${softLinkFolder}/${softLinkBAI}

# use relative path to make the softlink still be available even the folder been moved to different place
if ${bContainUMIReads}; then
  umiBAM=`basename ${OUTBAM_DEDUP_NOPHIX_UMI}`
  umiBAI=`basename ${OUTBAI_DEDUP_NOPHIX_UMI}`
  cd ${softLinkFolder} && ln -s ${umiBAM} ${softLinkBAM}
  cd ${softLinkFolder} && ln -s ${umiBAI} ${softLinkBAI}
else
  nophixBAM=`basename ${OUTBAM_NOPHIX}`
  nophixBAI=`basename ${OUTBAI_NOPHIX}`
  cd ${softLinkFolder} && ln -s ${nophixBAM} ${softLinkBAM}
  cd ${softLinkFolder} && ln -s ${nophixBAI} ${softLinkBAI}
fi

# make output files group RW-able
chmod g+rw ${OUTDIR}/${OUTNAME_SIMPLE}*.*
chmod g+rw $FLAGSTAT_FILE

if [[ ${iMergedSample} -eq 1 ]]; then
    rm -f ${FASTQ1}
    rm -f ${FASTQ2}
fi

# set flags
if [[ -f ${flagWorking} ]];  then
  rm ${flagWorking}
fi
touch ${flagDone}

echo ====================
date
echo Done!
echo ====================