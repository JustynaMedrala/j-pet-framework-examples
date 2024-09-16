#!/bin/bash
i=0

rm /mnt/home/jmedrala/output/cat_data/cosmics/*
for f in /mnt/jpet/Preselection/Run11/data/2020.04.09_0002/*.presel.evt.root; #cosmics
do 
	if [ "$i" = 500 ]; then break; fi
	sbatch run_L_cosmics.sh $f
	i=$((i+1))
done
