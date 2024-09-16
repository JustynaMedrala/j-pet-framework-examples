#!/bin/bash
#SBATCH --job -name=framework_MC
#SBATCH --partition=INTEL_HASWELL
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

if [ -z ${1+x} ]; then echo "give a file name as argument"; exit; fi

/mnt/home/jmedrala/j-pet-framework-examples/build/LargeBarrelAnalysis/LargeBarrelAnalysis.x -t root -f $1 -u /mnt/home/jmedrala/j-pet-framework-examples/LargeBarrelAnalysis/userParams_mc.json -i 9 -l /mnt/home/jmedrala/j-pet-framework-examples/CalibrationFiles/9_RUN/detectorSetupRun9.json -o /mnt/home/jmedrala/output/cat_mc/$2
