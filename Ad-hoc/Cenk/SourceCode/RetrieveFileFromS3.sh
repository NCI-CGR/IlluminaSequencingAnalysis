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
    fi
    argIndex=$((argIndex + 1))
done

SECONDS=0

# Get BAM
CMD="obj_get -v DCEG_COVID_WGS ${strBAMFileS3} -D ${strDestDir} -p -V --strip 12"
echo "${CMD}"
eval "${CMD}"

result=$?
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

if [[ ${result} -eq 0 ]]; then
  echo "$(date): ${strDestFileName} has been reteieved successfully!"
else
  echo "Error: $(date) ${strDestFileName} was failed to be retrieved!"
  exit 1
fi

echo 

# Get BAI
CMD="obj_get -v DCEG_COVID_WGS ${strBAMFileS3}.bai -D ${strDestDir} -p -V --strip 12"
echo "${CMD}"
eval "${CMD}"

# Check Running Results
result=$?
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

if [[ ${result} -eq 0 ]]; then
  echo "$(date): BAI has been reteieved successfully!"
else
  echo "Error: $(date) BAI was failed to be retrieved!"
  exit 1
fi