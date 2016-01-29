#!/bin/bash

#$ -t 1-2
#$ -cwd
#$ -N homophy
#$ -j y

genes=("SEQ1" "SEQ2")
gene=${genes[$(expr $SGE_TASK_ID - 1)]}

if [ -e data/seed/$gene.fa ]; then
	make GENE=$gene DB=0
fi

