
DATE=`date '+%Y-%m-%d-%H-%M'`


CMD="${SCRIPT_DIR}/misc_merge_all.sh ${BUFFER_DIR}/output/bam_results_${DATE} ${BUFFER_DIR}/output/bam_results_${DATE}/all_out.txt"
rm -f ${BUFFER_DIR}/output/bam_results_${DATE}/all_out.txt
echo $CMD

CMD="java -jar ${SCRIPT_DIR}/secondaryParsing.jar compare $MANIFEST $BAM_OUT $REPORT_OUT $REPORT_ERR"
echo $CMD
eval $CMD
