#!/bin/bash

db_dir=$1

seqs=(`ls $db_dir`)

for seq in  ${seqs[@]}; do
	formatdb -i $db_dir/$seq -p T -l formatdb.log
done
