#!/bin/bash 

cat all.108.139.247.tumor.rsem.genes.results.txt | tail -n +2 | awk '{print $1}' | awk -F '|' '{print $1}' | sort > ./all.sort.gene.txt

cat all.108.139.247.tumor.rsem.genes.results.txt | tail -n +2 | awk '{print $1}' | awk -F '|' '{print $1}' | sort | uniq > ./all.sort.uniq.gene.txt

diff ./all.sort.gene.txt ./all.sort.uniq.gene.txt
