#!/bin/bash
#SBATCH --partition=INTEL_HASWELL
#SBATCH --job -name=framework_MC
#SBATCH --time=05:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

#echo "/mnt/home/jmedrala/output/unk/$2_larger/"

/mnt/home/jmedrala/j-pet-framework-examples/build/MCGeantAnalysis/MCGeantAnalysis.x -t mcGeant -f $1 -u /mnt/home/jmedrala/j-pet-framework-examples/MCGeantAnalysis/userParams.json -i 9 -l /mnt/home/jmedrala/j-pet-framework-examples/CalibrationFiles/9_RUN/detectorSetupRun9.json -o /mnt/home/jmedrala/output/unk/

#rm $1
base=$(basename ${1})
#rm /mnt/home/jmedrala/output/unk/$2_larger/${base%.*.*}.hits.evt.root
