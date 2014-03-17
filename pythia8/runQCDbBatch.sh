#!/bin/bash
#
#PBS -l walltime=03:00:00
#
cd ${HOME}/pythia8183/examples
${HOME}/pythia8183/examples/mainQCDb.exe ${name}
