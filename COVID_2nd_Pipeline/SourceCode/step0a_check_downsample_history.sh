#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")
MANIFEST=$1
DIR=`dirname $MANIFEST`
BASE_NAME=`basename $MANIFEST .csv`
tr -d '\15\32' < $MANIFEST > ${DIR}/${BASE_NAME}_unix.csv
mkdir -p ${DIR}/../variant_calling_logs/downsample
while IFS=$',' read -r -a myArray
do
   # echo ${myArray[14]} 
  if [[ "${myArray[5]}" != "CGF"* ]]; then
    if [[ "${myArray[0]}" == "E0296"* ]]; then
      PLATFORM="NextSeq"

    else
      if [[ "${myArray[0]}" == "E0265"* ]]; then
        PLATFORM="MiSeq"
      else
        PLATFORM="HiSeq"
      fi
    fi
    
    SEQ_DIR="/DCEG/CGF/Sequencing/Illumina/${PLATFORM}/PostRun_Analysis/Data"
    BAM_FILE=${SEQ_DIR}/*${myArray[2]}/BAM/${myArray[5]}_${myArray[4]}_L00${myArray[3]}_HQ_paired_dedup_properly_paired_nophix.bam
    for NAME in ${SEQ_DIR}/*${myArray[2]}/BAM/${myArray[5]}_${myArray[4]}_L00${myArray[3]}*downsample*; do
        if [[ -f "$NAME" ]]; then
            echo "Warning: ${SEQ_DIR}/*${myArray[2]}/BAM/${myArray[5]}_${myArray[4]}_L00${myArray[3]}.bam is downsampled!"
        fi
    done
         # echo $BAM_FILE
         # rm -f ${DIR}/../variant_calling_logs/downsample/${myArray[2]}_${myArray[5]}_${myArray[4]}_L00${myArray[3]}_downsampling.std???


  fi

done < ${DIR}/${BASE_NAME}_unix.csv
