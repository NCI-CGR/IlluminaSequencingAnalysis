#!/bin/bash

set -o pipefail

argIndex=1
for arg in "$@"
do
    if [[ ${argIndex} == 1 ]]; then
        strBAMFileS3=${arg}
    elif [[ ${argIndex} == 2 ]]; then
        strDestFileName=${arg}    
    elif [[ ${argIndex} == 3 ]]; then
        strDestDir=${arg}
    elif [[ ${argIndex} == 4 ]]; then
        coreNum=${arg}        
    fi
    argIndex=$((argIndex + 1))
done

SECONDS=0

strBEDFile="../Files/igh_orphons_GRCh38.bed"

strDestBAM="${strDestDir}"/"${strDestFileName}"

# Run Samtools view
CMD="samtools view -@ ${coreNum} -h -L ${strBEDFile} -o ${strDestBAM} ${strBAMFileS3} && samtools index -@ ${coreNum} ${strDestBAM}"

echo ${CMD}
eval ${CMD}

echo 

# Check Running Results
result=$?
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

if [[ ${result} -eq 0 ]]; then
  echo "$(date): ${strDestFileName} has been processed (target BAM) successfully!"
else
  echo "Error: $(date) ${strDestFileName} was failed to be processed (target BAM)!"
  exit 1
fi