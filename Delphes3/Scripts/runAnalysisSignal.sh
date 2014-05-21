#!/bin/bash
#PBS -l walltime=6:00:00
cd ${HOME}/Delphes-3.0.12
touch examples/mainAnalysis.cpp
make
time mainAnalysis --source signal --swapMuRandomly true --doHLT true
