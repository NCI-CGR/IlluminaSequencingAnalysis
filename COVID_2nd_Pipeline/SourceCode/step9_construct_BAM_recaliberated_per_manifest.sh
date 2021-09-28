#!/bin/sh

SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")
. ${SCRIPT_DIR:-.}/global_config_bash.rc
if [[ $# -ne 3 ]]; then
   echo "Usage: step9_construct_BAM_recaliberated_per_manifest.sh $Manifest_File $Build_name SOMATIC|GERMLINE"
   echo "Example: step9_construct_BAM_recaliberated_per_manifest.sh  /DCEG/Projects/Exome/builds/build_SR0493-001_Data_Delivery_2019_22641/Manifest/SR0493-001-ANALYSIS-MANIFEST-4-4-2019.csv build_SR0493-001_Data_Delivery_2019_22641 SOMATIC"
   exit 1 
fi

MANIFEST=$1
BUILD_NAME=$2
TYPE=$3

if [[ ! -f $MANIFEST ]] ; then
  echo "Error: cannot find manifest file, please check the directory of manifest file."
  exit 1
fi
BAM_BUILD=/DCEG/Projects/Exome/builds/${BUILD_NAME}/bam_location
BAM_MISSING_LIST=/DCEG/Projects/Exome/builds/${BUILD_NAME}/missing_bam_file.txt
rm -f $BAM_MISSING_LIST
mkdir -p $BAM_BUILD
total_samples=0
count_for_missing=0
count_for_present=0

for ID in $(awk -F "," '{if (NR > 1) {print $7"_"$13}}' $MANIFEST | sort | uniq); do
#or ID in $(awk '{print $1}' $Manifest | sort | uniq); do 
	GROUP=$(echo $ID | cut -d_ -f1)
	echo $ID
	echo $GROUP
        if [[ $TYPE == "SOMATIC" ]]; then
           BAM_INPUT=${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${ID}_bqsr_final_out.bam  
           BAI_INPUT=${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${ID}_bqsr_final_out.bam.bai     
        else
	   BAM_INPUT=${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${ID}.bam
	   BAI_INPUT=${BAM_REFORMATTED_RECALIBRATED_DIR}/${GROUP}/${ID}.bam.bai
        fi
        BAM_OUTPUT=${BAM_BUILD}/${GROUP}/${ID}.bam
        BAI_OUTPUT=${BAM_BUILD}/${GROUP}/${ID}.bam.bai
	echo $BAM_INPUT
	echo $BAM_OUTPUT
	#echo $BAI_INPUT	

	if [[ ! -f $BAM_INPUT ]]; then
		echo "Error: the BAM file $BAM_INPUT is missing!"
		let "count_for_missing=count_for_missing+1"
		echo $BAM_INPUT >> $BAM_MISSING_LIST
		#exit 1
	else
                if [[ ! -d ${BAM_BUILD}/${GROUP} ]]; then
		  CMD="mkdir -p ${BAM_BUILD}/${GROUP}"
                  echo $CMD
	    	  eval $CMD
                fi
		CMD="ln -s $BAM_INPUT $BAM_OUTPUT; ln -s $BAI_INPUT $BAI_OUTPUT"
                echo $CMD
		eval $CMD
		let "count_for_present=count_for_present+1"
	fi
	let "total_samples=total_samples+1"
	echo
done

chmod -R 775 $BAM_BUILD
echo
echo
echo "Total unique samples inside manifest file = $total_samples"
echo "Total missing samples from Bam recaliberated directory = $count_for_missing"
echo "Total samples linked from Bam recaliberated directory = $count_for_present"
echo
echo
