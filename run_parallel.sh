#!/bin/bash

#$ -t 1-2
#$ -cwd
#$ -N homophy
#$ -j y

label=$1
genes=("SEQ1" "SEQ2")
gene=${genes[$(expr $SGE_TASK_ID - 1)]}

if [ -e data/$label/$gene.fa ]; then
	make GENE=$gene DB=0 LAB=$label HO=yes
fi

