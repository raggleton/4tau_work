#!/bin/bash
#
#PBS -l walltime=11:00:00
#
cd ${HOME}/4Tau/CALCHEP/
${HOME}/4Tau/CALCHEP/calchep_batch ${filename}
