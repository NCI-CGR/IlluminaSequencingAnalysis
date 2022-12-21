#!/bin/bash

. /etc/profile.d/modules.sh
module purge

module load miniconda/4.8.3 sge

#. /home/lix33/.bashrc

#module load sge
#unset module
#source /scratch/lix33/lxwg/software/conda/miniconda3/condabin/activate envCGR

source /home/lix33/.conda/envs/envCGR/etc/profile.d/conda.sh
conda activate /home/lix33/.conda/envs/envCGR

#module load python3/3.8.12 samtools/1.8 tabix/1.9 bgzip/0.2.6 gcc/7.2.0 zlib/1.2.11 bedtools/2.27.1 python3/3.7.0 
module load samtools/1.8 tabix/1.9 bgzip/0.2.6 gcc/7.2.0 zlib/1.2.11 bedtools/2.27.1 python3/3.7.0 python3/3.8.12

ref=$1
reads=$2
sampleName=$3
outputBAMDir=$4
outputFlagDir=$5
threads=$6

echo "ref      : ${ref}"
echo "reads    : ${reads}"
echo "outputDir: ${outputBAMDir}"

flagWorking="${outputFlagDir}/${sampleName}.mapping.working"
flagDone="${outputFlagDir}/${sampleName}.mapping.done"

touch ${flagWorking}

outputBAMRegDir="${outputBAMDir}/Raw"
if [[ ! -d ${outputBAMRegDir} ]]; then
    mkdir -p ${outputBAMRegDir}
fi

OutSortBAM="${outputBAMRegDir}/${sampleName}.sorted.bam"
SortedIndex="${OutSortBAM}.bai"

if [[ ! -f ${OutSortBAM} ]] || [[ ! -f ${SortedIndex} ]]; then
    # Alignment
    OutSAM="${outputBAMRegDir}/${sampleName}.sam"
    CMD="minimap2 -t ${threads} -Y --MD -a ${ref} ${reads} > ${OutSAM}"
    echo ${CMD}
    echo
    eval ${CMD}
    if [[ $? -ne 0 ]]; then
        echo "Error: $(date) minimap2 is failed!"
        exit 1
    fi
    
    # SAM to BAM
    OutBAM="${outputBAMRegDir}/${sampleName}.bam"
    CMD="samtools view -@ ${threads} -S -b ${OutSAM} > ${OutBAM}"
    echo ${CMD}
    echo
    eval ${CMD}
    if [[ $? -ne 0 ]]; then
        echo "Error: $(date) samtools view is failed!"
        exit 1
    fi
    
    # Sort 
    OutSortBAM="${outputBAMRegDir}/${sampleName}.sorted.bam"
    CMD="samtools sort -@ ${threads} -o ${OutSortBAM} ${OutBAM}"
    echo ${CMD}
    echo
    eval ${CMD}
    if [[ $? -ne 0 ]]; then
        echo "Error: $(date) samtools sort is failed!"
        exit 1
    fi
    
    # Build Index
    SortedIndex="${OutSortBAM}.bai"
    CMD="samtools index -@ ${threads} ${OutSortBAM} ${SortedIndex}"
    echo ${CMD}
    echo
    eval ${CMD}
    if [[ $? -ne 0 ]]; then
        echo "Error: $(date) samtools index is failed!"
        exit 1
    fi
    
    # Remove SAM, BAM
    rm ${OutSAM}
    rm ${OutBAM}
fi

# Extract Mapped Reads
outputBAMMappedDir="${outputBAMDir}/Mapped"
if [[ ! -d ${outputBAMMappedDir} ]]; then
    mkdir -p ${outputBAMMappedDir}
fi
    
mappedReadsBAM="${outputBAMMappedDir}/${sampleName}.mapped.sorted.bam"
CMD="samtools view -@ ${threads} -q 10 -F 4 -F 256 -h ${OutSortBAM} | grep -v -E -e '\bXA:Z:' -e '\bSA:Z:' | samtools view -b - > ${mappedReadsBAM}"
echo
echo "Extract Mapped BAM -->"
echo ${CMD}
eval ${CMD}

SortedIndex="${mappedReadsBAM}.bai"
CMD="samtools index -@ ${threads} ${mappedReadsBAM} ${SortedIndex}"
echo
echo ${CMD}
eval ${CMD}

# Get gene list
listGene="${mappedReadsBAM}.gene.list"
CMD="samtools view ${mappedReadsBAM} | awk '{print $1}' > ${listGene}"
echo
echo ${CMD}
eval ${CMD}

# Extract Unmapped reads
outputBAMUnMappedDir="${outputBAMDir}/UnMapped"
if [[ ! -d ${outputBAMUnMappedDir} ]]; then
    mkdir -p ${outputBAMUnMappedDir}
fi

unmappedReadsBAM="${outputBAMUnMappedDir}/${sampleName}.unmapped.bam"
CMD="samtools view -@ ${threads} -h -b -f 4 ${OutSortBAM} > ${unmappedReadsBAM}"
echo
echo "Extract UnMapped BAM -->"
echo ${CMD}
eval ${CMD}

bamIndex="${unmappedReadsBAM}.bai"
CMD="samtools index -@ ${threads} ${unmappedReadsBAM} ${bamIndex}"
echo
echo ${CMD}
eval ${CMD}

# Get gene list
listGene="${unmappedReadsBAM}.gene.list"
CMD="samtools view ${unmappedReadsBAM} | awk '{print $1}' > ${listGene}"
echo
echo ${CMD}
eval ${CMD}

rm ${flagWorking}
touch ${flagDone}

