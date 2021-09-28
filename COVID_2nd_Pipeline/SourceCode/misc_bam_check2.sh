#!/bin/sh
BAM_FILE=$1
OUTPUT_FILE=$2
BASE_NAME=$3
OUT_DIR=`dirname $OUTPUT_FILE`

set -o pipefail
module load samtools
TMP_DIR=/ttemp/wmy

# CMD="samtools view $BAM_FILE | awk -F\":\" '{printf(\"%s\t%s\t%s\t%s\n\",\$1,\$3,\$4,\$NF);}' | sort -u > $OUTPUT_FILE"
CMD="samtools view $BAM_FILE | awk -F\"\t\" '{for (i=1;i<=NF;i++) {if (\$i~/^RG/) {split(\$1,aa,\":\"); print aa[1]\"\t\"aa[3]\"\t\"aa[4]\"\t\"$i}}}' | sort -u > \${TMP_DIR}/\${BASE_NAME}_list"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
  echo "Error: samtools for reads is failed!"
  exit 1
fi
CMD="samtools view -H $BAM_FILE | awk -F\"\t\" '\$1~/@RG/ {for (i=1;i<=NF;i++) {if (\$i~/^ID:/) split(\$i,aa,\":\"); if (\$i~/^PU:/) split(\$i,bb,\".\");} print aa[2]\"\t\"bb[3]}' >  \${TMP_DIR}/\${BASE_NAME}_header"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
  echo "Error: samtools for the header is failed!"
  exit 1
fi


C1=`cat ${TMP_DIR}/${BASE_NAME}_header | wc -l`
C2=`cat ${TMP_DIR}/${BASE_NAME}_list | wc -l`
if [[ $C1 -ne $C2 ]]; then
  echo "Error: lane numbers are not equal!"
  exit 1
fi

awk -F"\t" '{
  if (FNR==NR){
  aa[$1]=$2
}
else{
  split($4,bb,":");
  print $1"\t"$2"\t"$3"\t"aa[bb[3]];
}
}' ${TMP_DIR}/${BASE_NAME}_header ${TMP_DIR}/${BASE_NAME}_list > $OUTPUT_FILE

if [[ $? -ne 0 ]]; then
  echo "Error: failed for samtools view!"
fi
rm -f ${TMP_DIR}/${BASE_NAME}_header ${TMP_DIR}/${BASE_NAME}_list
touch ${OUT_DIR}/${BASE_NAME}.DONE
