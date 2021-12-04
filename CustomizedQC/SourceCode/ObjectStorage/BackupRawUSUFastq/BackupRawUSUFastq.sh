#!/bin/bash

set -o pipefail

strRootDir="/data/COVID_WGS/COVNET_Data_Delivery/NCI_COVNET_WGS/Archive"

CMD="find ${strRootDir} -iname '*.fastq.gz' | \
     while read line; do obj_put -v DCEG_COVID_WGS -p \"lix33/COVID19/USU/RAW/Fastq\" ${line} -V -F --dry-run; done"
                                    
echo ${CMD}
echo "---"

SECONDS=0

echo -e "Run Command line \n >>>>>"

find ${strRootDir} -iname '*.fastq.gz' | while read line; do obj_put -v DCEG_COVID_WGS -p "lix33/COVID19/USU/RAW/Fastq" ${line} -V -F; done

echo -e "<<<<<< \n End!" 
       
result=$?
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"

if [[ ${result} -eq 0 ]]; then
  echo "$(date): ${strUSUSampleID} in ${strFlowcellName} has been reteieved successfully!"
  
  rm ${strFlagWorking}
  touch ${strFlagDone}
else
  echo "Error: $(date) ${strUSUSampleID} in ${strFlowcellName} was failed to be retrieved!"
  exit 1
fi