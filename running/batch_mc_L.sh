#!/bin/bash

count=0
for f in /mnt/home/jmedrala/output/unk/$1_larger/*.unk.evt.root; 
#for f in /mnt/home/ptanty/tmp/*.unk.evt.root;
do
	base=$(basename ${f})
        file_name=${base%.*.*.*}
        #if [ "$count" = 1 ]; then break; fi
        if [ "$(ls -l /mnt/home/jmedrala/output/cat_mc/$1/${file_name}.cat.evt.root | wc -l)" == "0" ]; then
	sbatch run_L.sh $f $1 
	let count++
	fi
done

#sleep 200
#rm /mnt/home/jmedrala/output/cat_mc/$1_larger/*.presel.evt.root
