#!/bin/bash

for f in /mnt/home/jmedrala/output/*.unk.evt.root; 
do 
	echo "Processing $f file..."; 
        sbatch run_F.sh $f
done
