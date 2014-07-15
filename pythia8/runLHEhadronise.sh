#!/bin/bash
#
#PBS -l walltime=6:59:00
#
cd ${HOME}/pythia8183/examples
${HOME}/pythia8183/examples/mainLHEhadronise.exe ${LHEname} ${hepmcname}
