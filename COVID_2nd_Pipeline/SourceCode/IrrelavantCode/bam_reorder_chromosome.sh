#!/bin/sh

# Change the order of chromosomes in legacy BAM files to GATK-bundle-compatible
# The chromosome order in the GATK-bundle is: M, 1...9, 10...19, 20, 21, 22, X, Y
SCRIPT=$(readlink -f "$0")
DCEG_SEQ_POOL_SCRIPT_DIR=$(dirname "$SCRIPT")
. ${DCEG_SEQ_POOL_SCRIPT_DIR:-.}/global_config_bash.rc

set -o pipefail
INPUT=$1
OUTPUT=$2

if [[ "$OUTPUT" == "" ]]; then
	OUTPUT=$INPUT
fi

OUTDIR=`dirname $OUTPUT`
OUTNAME=`basename $OUTPUT .bam`
OUTBAI=$OUTDIR/${OUTNAME}.bam.bai

TMP_OUTPUT_NAME=$TMP_DIR/${OUTNAME}_reorder_chr

echo ==========
date
echo Change the listing order of chromosomes in the BAM file for compatibility to GATK bundle
echo Input: $INPUT
echo Output: $OUTPUT
echo ==========

echo "[$(date +%H:%M:%S)] Changing BAM chromosome listing order..."
(
# print the non-SQ line of the BAM header
samtools view -H $INPUT | awk '$1 != "@SQ"'
# get the @SQ lines of header, which are the chromosome information, 
# then extract the particular line according to the custom order below
for n in M $(seq 1 22) X Y; do 
#for n in $(cat ~/b37_ref/human_g1k_v37.dict | cut -f2 | cut -f2 -d: | tail -n +2); do
	samtools view -H $INPUT | awk -v n=$n '$1 == "@SQ" && $2 == "SN:chr"n' | sed s/16569/16571/g
	#samtools view -H $INPUT | awk -v n=$n '$1 == "@SQ" && $2 == "SN:"n'
done
# output the contents; making all information output so far a whole SAM stream but with changed chromosome order in the header
# this SAM stream is further changed to binary BAM format and re-sorted
samtools view $INPUT
) | samtools view -bS - | samtools sort -m 2000000000 -o $TMP_OUTPUT_NAME
# Patch the samtools stderr to indicate this portion of work is done.

if [[ $? -ne 0 ]]; then
   echo "Error: Changing BAM chromosome order is failed!"
   exit 1
else
   echo "Done." 1>&2
fi
# move the ordered file to final destination
# remove old indices, and create new index
echo "[$(date +%H:%M:%S)] Doing some more housekeeping work..."

CMD="mv ${TMP_OUTPUT_NAME}.bam $OUTPUT"
echo "[$(date +%H:%M:%S)] Command: $CMD"
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: $CMD is failed!"
   exit 1
fi
CMD="rm -rf $OUTDIR/$OUTNAME.bam.bai $OUTDIR/$OUTNAME.bai; samtools index $OUTPUT $OUTBAI"
echo "[$(date +%H:%M:%S)] Command: $CMD"
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: $CMD is failed!"
   exit 1
fi
# patching the permissions of the files
chgrp ncicgf_dceg_exome $OUTPUT $OUTBAI 2>/dev/null
chmod og+r-x $OUTPUT $OUTBAI 2>/dev/null

echo ==========
date
echo Done!
echo ==========
