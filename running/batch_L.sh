#!/bin/bash
i=0

#rm /mnt/home/jmedrala/output/cat_data/cosmics/*.cat*
#for f in /mnt/jpet/Preselection/Run11/data/2020.04.09_0002/*.presel.evt.root; #cosmics
#rm /mnt/home/jmedrala/output/cat_data/*.cat*
for f in /mnt/jpet/Preselection/Run11/data/2020.04.16_0001/*.presel.evt.root; 
#for f in /mnt/jpet/Preselection/Run11/data/2020.12.30_0001/*.presel.evt.root; 
do 
	if [ "$i" == 4000 ]; then break; fi
        sbatch run_L.sh $f
	i=$((i+1))
done
