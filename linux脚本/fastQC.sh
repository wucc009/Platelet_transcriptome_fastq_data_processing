#!/bin/bash
# ------------------------------------
#              fastQC
# ------------------------------------
echo -e "\n\n\n\e[32m\e[1m------fastQC------\e[0m\n\n\n"
mkdir -p ./result/fastQC/$1
fastqc -q -t 5 -o ./result/fastQC/$1 ./rawdata/FTP/$1/*gz # -q静默状态
# ------------------------------------
#               MultiQC
# ------------------------------------
echo -e "\e[32m\e[1m------MultiQC------\e[0m\n\n\n"
mkdir -p ./result/MultiQC/$1
multiqc ./result/fastQC/$1/*.zip -o ./result/MultiQC/$1
