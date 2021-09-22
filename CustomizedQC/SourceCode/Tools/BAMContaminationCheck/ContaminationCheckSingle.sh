#/bin/bash

argIndex=1
for arg in "$@"
do
    if [[ ${argIndex} == 1 ]]; then
        vcfFile=${arg}
    elif [[ ${argIndex} == 2 ]]; then
        bamFile=${arg}
    elif [[ ${argIndex} == 3 ]]; then
        outDir=${arg}
    elif [[ ${argIndex} == 4 ]]; then
        flagWorking=${arg}
    elif [[ ${argIndex} == 5 ]]; then
        flagDone=${arg}
    fi
    argIndex=$((argIndex + 1))
done

#load module
module load verifybamid/1.1.3
which verifyBamID 


CMD="verifyBamID --ignoreRG --chip-none --precise --best --vcf ${vcfFile} --bam ${bamFile} --out ${outDir}"
echo "Check BAM contamination!"
echo 

echo ${CMD}
echo

SECONDS=0
eval ${CMD}
result=$?
duration=$SECONDS
echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
if [[ ${result} -ne 0 ]]; then
  echo "Error: $(date) BAM contamination check was failed!"
  exit 1
else
  echo "$(date) BAM contamination check was finished successfully!"
fi

#Change flag
if [ -f ${flagWorking} ]; then
    rm ${flagWorking}
fi

touch ${flagDone}
