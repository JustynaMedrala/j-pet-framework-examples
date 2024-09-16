#!/bin/bash

if [ -z ${1+x} ]; then echo "give a file name as argument"; exit; fi

/mnt/home/jmedrala/j-pet-framework-examples/build/FlatTree/FlatTree.x -t root -f $1 -u /mnt/home/jmedrala/j-pet-framework-examples/LargeBarrelAnalysis/userParams.json -i 9 -l /mnt/home/jmedrala/j-pet-examples/CalibrationFiles/9_RUN/detectorSetupRun9.json
