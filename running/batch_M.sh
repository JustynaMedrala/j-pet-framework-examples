#!/bin/bash
count=0
for f in /mnt/home/jmedrala/output/mcGeant/$1_larger/*.mcGeant.root; 
do 
	base=$(basename ${f})
	file_name=${base%.*.*}
	if [ "$(ls -l /mnt/home/jmedrala/output/unk/$1_larger/${file_name}.unk.evt.root | wc -l)" == "0" ]; then
	#echo "Processing $f file..."; 
        sbatch run_M.sh $f $1
	let count++
	fi
        
done
echo "$count"
