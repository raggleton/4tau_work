#!/bin/bash
#
#PBS -l walltime=119:59:59
#
cd /panfs/panasas01/phys/ra12451/pythia8185/examples
mainQCD.exe ${args}
