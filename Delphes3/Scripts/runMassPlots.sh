#!/bin/bash
#PBS -l walltime=3:00:00
cd ${HOME}/Delphes-3.0.12
touch examples/massPlots.cpp
make
time massPlots --source qcdb --swapMuRandomly true --doHLT true 
