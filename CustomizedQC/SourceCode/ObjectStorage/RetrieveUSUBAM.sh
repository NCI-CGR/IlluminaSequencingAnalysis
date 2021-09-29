#!/bin/bash

argIndex=1
for arg in "$@"
do
    if [[ ${argIndex} == 1 ]]; then
        strFlowcellName=${arg}
    elif [[ ${argIndex} == 2 ]]; then
        strUSUSampleID=${arg}    
    elif [[ ${argIndex} == 3 ]]; then
        strDir=${arg}
    elif [[ ${argIndex} == 4 ]]; then
        strFlagWorking=${arg}
    elif [[ ${argIndex} == 5 ]]; then
        strFlagDone=${arg}
    fi
    argIndex=$((argIndex + 1))
done


CMD="obj_ls -v DCEG_COVID_WGS \
            -h \
            -m *${strFlowcellName}*${strUSUSampleID}*dedup_nophix.???  
            |awk 'NR>1 {print \$8}' |  
            while read line; do obj_get -v DCEG_COVID_WGS ${line} 
                                        -D ${strDir} 
                                        -p -V --strip 6 --dry-run; done"
                                    
echo ${CMD}
echo "---"
SECONDS=0

echo -e "Run Command line \n >>>>>"

obj_ls -v DCEG_COVID_WGS -h -m *${strFlowcellName}*${strUSUSampleID}*dedup_nophix.??? \
| awk 'NR>1 {print $8}' \
| while read line; do obj_get -v DCEG_COVID_WGS ${line} -D ${strDir} -p -V --strip 6; done

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