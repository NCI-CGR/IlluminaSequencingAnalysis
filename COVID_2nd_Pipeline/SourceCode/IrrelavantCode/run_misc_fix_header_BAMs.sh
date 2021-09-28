#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_HOME=$(dirname "$SCRIPT")

MANIFEST=/DCEG/Projects/Exome/builds/2016-09-07_with_population_control/manifest/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv
if [[ ! -f $MANIFEST ]] ; then
	echo "Error: cannot find manifest file, please check the directory of manifest file."
	exit 1
fi
BAM_DIR=/DCEG/Projects/Exome/SequencingData/BAM_original
OUT_DIR=/DCEG/Projects/Exome/SequencingData/BAM_processed
LOG_DIR=${OUT_DIR}/LOGs
mkdir -p $OUT_DIR
mkdir -p $LOG_DIR
for ID in $(awk -F "," '{if (NR > 1) {print $7"_"$13}}' $MANIFEST | sort | uniq); do
#or ID in $(awk '{print $1}' $Manifest | sort | uniq); do 
        
  GROUP=$(echo $ID | cut -d_ -f1)
  ANALYSIS_ID=`echo $ID | awk -F"_" '{{for (i=2;i<NF;i++) printf $i"_"} print $NF}'`
  # echo $ID"\t"$GROUP"\t"$ANALYSIS_ID 
  BAM_INPUT=${BAM_DIR}/${GROUP}/${ID}.bam
  BAM_OUTPUT=${OUT_DIR}/${GROUP}/${ID}.bam
  if [[ "$ANALYSIS_ID" == "WM_4954_1003_A" || "$ANALYSIS_ID" == "ACS_33267_CO_U" || $ANALYSIS_ID == "BC_3113_2002_A" || "$ANALYSIS_ID" == "BLC_BPLC_313864_A" ]]; then
    rm -f ${LOG_DIR}/fix_${ANALYSIS_ID}.std??? 
    CMD="qsub -q seq*.q,long.q -o ${LOG_DIR}/fix_${ANALYSIS_ID}.stdout -e ${LOG_DIR}/fix_${ANALYSIS_ID}.stderr -N ${ANALYSIS_ID}_fix -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/misc_fix_bam_header.sh $BAM_INPUT $BAM_OUTPUT $GROUP $ANALYSIS_ID $MANIFEST"
    echo $CMD
    eval $CMD
  fi

done
#                for LINE in $(awk -F "," -v analysis_id=$ID '{if (NR > 1) && ($13==analysis_id) {print NR}}' $Manifest); do
#                   echo $LINE
#                   RG_ID=
#                   RG_LB=
#                   RG_PL=
#                   RG_SM=
#                   RG_PU=
#		cmd="mkdir -p $BAM_recaliberated_per_manifest/$GROUP 2>/dev/null"
