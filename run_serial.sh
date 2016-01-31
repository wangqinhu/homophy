#!/bin/bash

label=$1
genes=("SEQ1" "SEQ2")

for gene in  ${genes[@]}; do
	if [ -e data/$label/$gene.fa ]; then
		make GENE=$gene DB=0 LAB=$label
	fi
done

./homostat.pl data/$label\_conf/$label\1.conf data/$label\_family > $label.tsv
