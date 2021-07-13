#!/bin/bash

source /etc/profile.d/modules.sh
module load python/3.7

python3 /home/lix33/lxwg/Git/IlluminaSequencingAnalysis/CustomizedQC/SourceCode/CustomizedQC.py /data/COVID_WGS/primary_analysis/COVID19/06_30_2021/ProcessedData
