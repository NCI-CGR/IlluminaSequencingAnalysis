#/bin/bash

SCRIPT="/home/lixin/lxwg/ad-hoc/Ecotyper/ecotyper/EcoTyper_recovery_bulk.R"
GENEExpression="/home/lixin/lxwg/ad-hoc/Ecotyper/NewData/Organized/HKBCBatch123/all.108.139.247.tumor.rsem.genes.results.no.dup.txt"
DIROutput="/home/lixin/lxwg/ad-hoc/Ecotyper/NewData/Organized/Run/EcoTyper_recovery_bulk/ORG_Discovery_Bulk_RNA"
ANNOTATION="/home/lixin/lxwg/ad-hoc/Ecotyper/NewData/Organized/HKBCBatch123/all.108.139.247.tumor.rsem.genes.results.no.dup.annotation"
DISCOVERY="Carcinoma"
DIRWork="/home/lixin/lxwg/ad-hoc/Ecotyper/ecotyper"

cd ${DIRWork} && Rscript ${SCRIPT} -d ${DISCOVERY}\
	                           -m ${GENEExpression} \
				   -a ${ANNOTATION} \
				   -c Tissue \
				   -o ${DIROutput} \
				   -t 4

