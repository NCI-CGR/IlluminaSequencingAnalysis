#!/bin/sh
. ${SCRIPT_HOME:-.}/global_config_bash.rc
if [ $# -ne 1 ]; then
  DATE_NEWER="2012-01-01"
else
  DATE_NEWER=$1
  if [[ $DATE_NEWER =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
    echo $DATE_NEWER is a valid YYYY-MM-DD date
  else
    echo "Error: $DATE_NEWER is NOT a valid YYYY-MM-DD date!"
    exit 1
  fi
fi
rm -f ${BUFFER_DIR}/input/timestamp
touch ${BUFFER_DIR}/input/timestamp -d $DATE_NEWER
if [[ $? -ne 0 ]]; then
  echo "Error: timestamp is failed to crreated!"
  exit 1
fi
set -e
SEQ_DIRs=(/DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data /DCEG/CGF/Sequencing/Illumina/MiSeq/PostRun_Analysis/Data /DCEG/CGF/Sequencing/Illumina/NextSeq/PostRun_Analysis/Data /DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data /DCEG/CGF/Sequencing/Illumina/MiSeq/PostRun_Analysis/Data /DCEG/CGF/Sequencing/Illumina/NextSeq/PostRun_Analysis/Data)
PLATFORMs=(HiSeq MiSeq NextSeq HiSeq MiSeq NextSeq)
LOG_FILEs=(1*/*/_sequence_quality_trimming_*.stdout 1*/*/_sequence_quality_trimming_*.stdout  1*/*/_sequence_quality_trimming_*.stdout 2*/*/_sequence_quality_trimming_*.stdout 2*/*/_sequence_quality_trimming_*.stdout  2*/*/_sequence_quality_trimming_*.stdout)
DATE_OUTPUT=`date +%Y%m%d%H%M`
ERR_FILE=${BUFFER_DIR}/input/qt_errors_${DATE_OUTPUT}.txt
OUT_FILE=${BUFFER_DIR}/input/qt_all_${DATE_OUTPUT}.txt
set -e
rm -f $ERR_FILE $OUT_FILE
for ((jjj=0;jjj<${#SEQ_DIRs[@]};++jjj)); do
  CMD="cd ${SEQ_DIRs[$jjj]}"
  echo $CMD
  eval $CMD
  if [[ $? -ne 0 ]]; then
    echo "Error: ${SEQ_DIRs[$jjj]} does not existed!"
    exit 1
  fi
  
  echo "[$(date)]: finding the log files newer than $DATE_NEWER ..."
  FILE_COUNT=`find ${LOG_FILEs[$jjj]} -newer ${BUFFER_DIR}/input/timestamp | wc -l`
  if [[ $FILE_COUNT -eq 0 ]]; then
    echo "Warning: ${SEQ_DIRs[$jjj]} does not existed!"
    continue
  fi
  FILE_LIST=`find ${LOG_FILEs[$jjj]} -newer ${BUFFER_DIR}/input/timestamp`
  echo $CMD
  eval $CMD
  for i in $FILE_LIST; do
    TRIM_VER=0
    LOG_DIR=`dirname $i`
    if [[ $i == *"K00278"* ]]; then  # only check HiSeq 4000
      TRIM_VER=`grep trimmomatic-0.36.jar $i | wc -l` 
    fi
    # echo $TRIM_VER    
    if [[ $TRIM_VER -eq 0 &&  "$i" == *"K00278"* ]]; then
        echo "Hiseq 4000 -- old trimming: ${SEQ_DIRs[$jjj]}/${i}"
        FILEs=`awk '{if ($0~/Output File 1/)
                      output1=$4;
                  if ($0~/Output File 3/)
                      output3=$4;}
                  END{
                    printf "%s;%s", output1,output3;
                  }' $i`
        if [[ $LOG_DIR != *"build"* ]]; then
          if [[ -f ${LOG_DIR}/_CASAVA_demultiplexing.stderr ]]; then
              COUNT=`awk 'BEGIN {count=0} $0~/bcl2fastq v2.17.1.14/ {count=1} END {printf "%d",count;} ' ${LOG_DIR}/_CASAVA_demultiplexing.stderr`
              if [[ $COUNT -eq 0 ]];  then
                  echo "Error: Not in new bcl2fastq for Hiseq 4000: ${SEQ_DIRs[$jjj]}/${i}" >> $ERR_FILE
              else
                  FILE1=`echo $FILEs | cut -d";" -f1`
                  FILE2=`echo $FILEs | cut -d";" -f2`
                  if [[ ! -s $FILE1 ]] || [[ ! -s $FILE2 ]]; then
                      echo "Warning: Fastq files not existed or empty for Hiseq 4000: $FILE1 and $FILE2! Log file: ${SEQ_DIRs[$jjj]}/${i}" >> $ERR_FILE
                  fi
                  if [[ "$FILE1" == *"DCEG"* ]]; then
                      FLOWCELL=`echo $FILE1 | cut -d"/" -f9`
                      LANE=`echo $FILE1 | cut -d"/" -f11`
                      PROJECT=`echo $FILE1 | cut -d"/" -f12`
                      SAMPLE=`echo $FILE1 | cut -d"/" -f13`
                      #  because sample ID may contain _; so have to read from end using _ as delimiter
                      INDEX=`echo $FILE1 | cut -d"/" -f14 | rev | cut -d"_" -f6 | rev`
                  else
                      FLOWCELL=`echo $FILE1 | cut -d"/" -f8`
                      LANE=`echo $FILE1 | cut -d"/" -f10`
                      PROJECT=`echo $FILE1 | cut -d"/" -f11`
                      SAMPLE=`echo $FILE1 | cut -d"/" -f12`
                      INDEX=`echo $FILE1 | cut -d"/" -f13 | rev | cut -d"_" -f6 | rev`
                  fi
                  DATE=`date +%Y-%m-%d -r ${SEQ_DIRs[$jjj]}/${i}`
                  echo -e ${PLATFORMs[$jjj]}"\t"$FLOWCELL"\t"$INDEX"\t"$LANE"\t"$PROJECT"\t"$SAMPLE"\tOLD Trimmomatic\t" ${SEQ_DIRs[$jjj]}/${i}"\t"$DATE >> $OUT_FILE
              fi
          else
              echo "Error: No demultiplexing log file: ${SEQ_DIRs[$jjj]}/${i}" >> $ERR_FILE
          fi
        else
          echo "Error: Old Trimming in build folder: ${SEQ_DIRs[$jjj]}/${i}" >> $ERR_FILE         
        fi
    else
      if [[ $TRIM_VER -eq 0 ]]; then
          if [[ $jjj -lt 1 ]]; then
             echo "Hiseq 2500: ${SEQ_DIRs[$jjj]}/${i}"
          else
             if [[ $jjj -eq 1 ]]; then
                 echo "MiSeq : ${SEQ_DIRs[$jjj]}/${i}"
             else
                 echo "NextSeq: ${SEQ_DIRs[$jjj]}/${i}"
             fi
          fi
      else
          echo "HiSeq 4000 (new trimming): ${SEQ_DIRs[$jjj]}/${i}"
      fi

      if [[ ! -r ${SEQ_DIRs[$jjj]}/${i} ]]; then
         echo "Error: ${SEQ_DIRs[$jjj]}/${i} has no read permission!"
         continue 
      fi
      COUNT=`awk 'BEGIN {
                    count=0; 
                    count1=0; 
                    count2=0; 
                    line=0; 
                  } 
                  {
                  if ($0~/trimmomatic-0.36.jar/) {
                    count++; 
                    line=NR
                    if ($0~/TruSeq3-PE-2.fa/) 
                       count1++
                  }
                  else{
                    if (line>0){
                       if (NR==line+1)
                         if ($0~/successfully/)
                            count2++; 
                    }
                  }
                  if ($0~/Output File 1/)
                      output1=$4;
                  if ($0~/Output File 3/)
                      output3=$4;
                }
           END {
              if (count>0)
                  printf "%d_%d_%d:%s;%s;",count,count1,count2,output1,output3;
              else
                  printf "NA"
           }' $i`
      echo ${SEQ_DIRs[$jjj]}/${i}":$COUNT"
      if [[ $COUNT != "NA" ]]; then
          C1=`echo $COUNT | cut -d"_" -f1`   
          C2=`echo $COUNT | cut -d"_" -f2`
          C3=`echo $COUNT | cut -d"_" -f3 | cut -d":" -f1`
          FILEs=`echo $COUNT | cut -d":" -f2`
          if [[ $C1 -gt 1 ]]; then
              echo "Warning: Multiple trimmomatic: ${SEQ_DIRs[$jjj]}/${i} $C1" >> $ERR_FILE
          fi
          if [ $C1 -ne $C2 ] || [ $C1 -ne $C3 ] || [ $C2 -ne $C3 ]; then
              echo "Error: Trimmomatic and adapter sequence and successfully are unmatched: ${SEQ_DIRs[$jjj]}/${i} $C1 $C2 $C3" >> $ERR_FILE
          fi
          if [ $C1 -eq $C2 ] && [ $C1 -eq $C3 ] && [ $C1 -gt 0 ]; then
              FILE1=`echo $FILEs | cut -d";" -f1`
              FILE2=`echo $FILEs | cut -d";" -f2`
              if [[ ! -s $FILE1 ]] || [[ ! -s $FILE2 ]]; then
                  echo -e "Warning: Fastq files are not existed or empty for ${PLATFORMs[$jjj]}: $FILE1 and $FILE2! Log file: ${SEQ_DIRs[$jjj]}/${i}" >> $ERR_FILE
              fi
              if [[ "$FILE1" == *"DCEG"* ]]; then
                 FIELD_NUM=`echo $FILE1 | awk -F"/" '{print NF}'`
                 FLOWCELL=`echo $FILE1 | cut -d"/" -f9`
                 if [[ $FIELD_NUM -eq 13 ]]; then
                    PROJECT=`echo $FILE1 | cut -d"/" -f11`
                    SAMPLE=`echo $FILE1 | cut -d"/" -f12`
                    INDEX=`echo $FILE1 | cut -d"/" -f13 | rev | cut -d"_" -f6 | rev`
                    FULL_LANE=`echo $FILE1 | cut -d"/" -f13 | cut -d"_" -f3`
                    LANE=`echo ${FULL_LANE/00/}`
                 else
                    LANE=`echo $FILE1 | cut -d"/" -f11`
                    PROJECT=`echo $FILE1 | cut -d"/" -f12`
                    SAMPLE=`echo $FILE1 | cut -d"/" -f13`
                    INDEX=`echo $FILE1 | cut -d"/" -f14 | rev | cut -d"_" -f6 | rev`
                 fi
              else
                FIELD_NUM=`echo $FILE1 | awk -F"/" '{print NF}'`
                if [[ $FIELD_NUM -eq 12 ]]; then
                  FLOWCELL=`echo $FILE1 | cut -d"/" -f8`
                  FULL_LANE=`echo $FILE1 | cut -d"/" -f12 | cut -d"_" -f3`
                  PROJECT=`echo $FILE1 | cut -d"/" -f10`
                  SAMPLE=`echo $FILE1 | cut -d"/" -f11`
                  INDEX=`echo $FILE1 | cut -d"/" -f12 | rev | cut -d"_" -f6 | rev`
                  LANE=`echo ${FULL_LANE/00/}`
                else
                  FLOWCELL=`echo $FILE1 | cut -d"/" -f8`
                  LANE=`echo $FILE1 | cut -d"/" -f10`
                  PROJECT=`echo $FILE1 | cut -d"/" -f11`
                  SAMPLE=`echo $FILE1 | cut -d"/" -f12`
                  INDEX=`echo $FILE1 | cut -d"/" -f13 | rev | cut -d"_" -f6 | rev`
                fi
              fi
              DATE=`date +%Y-%m-%d -r ${SEQ_DIRs[$jjj]}/${i}`
              echo -e ${PLATFORMs[$jjj]}"\t"$FLOWCELL"\t"$INDEX"\t"$LANE"\t"$PROJECT"\t"$SAMPLE"\tNEW Trimmomatic\t"${SEQ_DIRs[$jjj]}/${i}"\t"$DATE >> $OUT_FILE 
          fi
      fi
    fi
  done
done
echo "Merging all old scan results ..."
CMD="cat ${BUFFER_DIR}/input/qt_all_base.txt ${BUFFER_DIR}/input/qt_all_${DATE_OUTPUT}.txt | sort | uniq > ${BUFFER_DIR}/input/qt_all.txt"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: merging is failed!"
   exit 1
fi
CMD="cp ${BUFFER_DIR}/input/qt_all.txt ${BUFFER_DIR}/input/qt_all_base.txt"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: cp is failed!"
   exit 1
fi
echo "[$(date)]: all done!"

