#!/bin/sh

# java -jar secondaryParsing.jar scan /DCEG/Projects/Exome/SequencingData/variant_scripts/logs/GATK /DCEG/Projects/Exome/SequencingData/variant_scripts/logs/all_merging_records.txt /DCEG/Projects/Exome/SequencingData/variant_scripts/logs/all_merging_err.txt
echo "Script Dir is: $SCRIPT_DIR"
. ${SCRIPT_DIR}/global_config_bash.rc

CMD="java -jar ${SCRIPT_DIR}/secondaryParsing.jar scan /DCEG/Projects/Exome/SequencingData/variant_scripts/logs/GATK ${BUFFER_DIR}/input/all_merging_records.txt ${BUFFER_DIR}/input/all_merging_err.txt"
echo $CMD
eval $CMD
#for Osteosarcoma build, the merge logs in 2014 has _P in the analysis ID, which caused the mismatch of the merge records for the following Osteosarcoma build (analysis ID does not has _P), remove _P in the logs to keep the consistency of analysis ID for the secondary analysis.
