#!/bin/bash
#
#PBS -l walltime=11:59:00
#
cd ${HOME}/pythia8185/examples
${HOME}/pythia8185/examples/mainLHEhadronise.exe ${LHEname} ${hepmcname}
