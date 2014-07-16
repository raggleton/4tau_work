#!/bin/bash
#
#PBS -l walltime=119:59:59
#
cd ${HOME}/pythia8185/examples
mainQCD.exe ${args}
