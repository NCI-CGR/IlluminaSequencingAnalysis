#!/bin/bash

YAML="/home/lixin/lxwg/ad-hoc/Ecotyper/NewData/Organized/Run/EcoTyper_discovery_bulk/config_discovery_bulk.yml"
SCRIPT="/home/lixin/lxwg/ad-hoc/Ecotyper/ecotyper/EcoTyper_discovery_bulk.R"
ROOTDIR="/home/lixin/lxwg/ad-hoc/Ecotyper/ecotyper"

cd ${ROOTDIR} && Rscript ${SCRIPT} -c ${YAML}
