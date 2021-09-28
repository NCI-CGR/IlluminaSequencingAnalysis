#!/bin/sh
BAM_FILE=$1
OUTPUT_FILE=$2
BASE_NAME=$3
set -o pipefail
set -e
module load samtools
OUT_DIR=`dirname $OUTPUT_FILE`
samtools view $BAM_FILE | awk -F":" '
               BEGIN{i=0;}
               {if(NR==1) 
                     {
                       MACHINE[i]=$1; 
                       FLOWCELL[i]=$3; 
                       LANE[i]=$4;
                       PU[i]=$NF}
                else  {
                    found=0;
                    for (j=0;j<=i;j++){
                      if (($1==MACHINE[j]) && ($3==FLOWCEL[j]) && (LANE[j]==$4) && ($NF==PU[j])) {
                        found=1;
                        break;  
                      }
                    }
                    if (found==0) {i++; MACHINE[i]=$1; FLOWCELL[i]=$3; LANE[i]=$4; PU[i]=$NF;} 
                }
                }
                END {for (j=0;j<=i;j++)
                      printf("%s\t%s\t%s\t%s\n",MACHINE[j],FLOWCELL[j],LANE[j],PU[j]);
                }'>$OUTPUT_FILE
touch ${OUT_DIR}/${BASE_NAME}.done
