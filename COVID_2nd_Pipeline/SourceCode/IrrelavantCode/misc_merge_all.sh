#!/bin/sh
DIR=$1
cd $DIR
if [[ ! -f  ${DIR}/all_out.txt ]]; then
  rm -f ${DIR}/all_out.txt
fi
for i in *.out; do
   echo $i
   if [[  -s $i ]]; then
      if [[ $i == *"_bqsr_final_out.out" ]]; then
        BASE_NAME=`basename $i _bqsr_final_out.out`
      else
        BASE_NAME=`basename $i .out`
      fi
      CMD="awk -v b=$BASE_NAME '{printf(\"%s\t%s\n\",b,\$0)};' $i >> ${DIR}/all_out.txt"
      echo $CMD
      eval $CMD
      if [[ $? -ne 0 ]]; then
        echo "Error: Merging is failed!"
        exit 1
      fi
   fi

done


CMD="java -jar ${SCRIPT_DIR}/secondaryParsing.jar compare ${DIR}/manifest.csv ${DIR}/all_out.txt ${DIR}/report_cmp.txt ${DIR}/err_cmp.txt"
echo $CMD
eval $CMD
