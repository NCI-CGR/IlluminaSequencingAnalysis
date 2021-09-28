#!/bin/sh


SCRIPT=$(readlink -f "$0")
SCRIPT_HOME=$(dirname "$SCRIPT")

. ${SCRIPT_HOME}/global_config_bash.rc
sed -i.bak 's/OSTE_P/OSTE/g' ${BUFFER_DIR}/input/all_merging_records.txt
