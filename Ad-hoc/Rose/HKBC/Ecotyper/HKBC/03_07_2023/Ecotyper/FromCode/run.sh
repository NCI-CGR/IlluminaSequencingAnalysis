#/bin/bash

SCRIPT="/home/lixin/lxwg/ad-hoc/Ecotyper/ecotyper/EcoTyper_recovery_bulk.R"
GENEExpression="/home/lixin/lxwg/ad-hoc/Ecotyper/HKBC/all.231tumor.normal.renamed.no.dup.txt"
DIROutput="/home/lixin/lxwg/ad-hoc/Ecotyper/HKBC/03_07_2023/Ecotyper/FromCode/Regular_Bulk_RNA"
ANNOTATION="/home/lixin/lxwg/ad-hoc/Ecotyper/HKBC/all.annotation.type.txt"

DIRWork="/home/lixin/lxwg/ad-hoc/Ecotyper/ecotyper"

cd ${DIRWork} && Rscript ${SCRIPT} -d Carcinoma -m ${GENEExpression} -a ${ANNOTATION} -c Tissue -o ${DIROutput} -t 4
