#!/bin/bash

module load python/3.7
module load samtools/1.13

keytable=$1
echo "keytable: "${keytable}

ABSOLUTE_PATH=$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")")
echo ${ABSOLUTE_PATH}

script="${ABSOLUTE_PATH}/AutoFramework.py"

python ${script} ${keytable}