#!/bin/bash
#
#PBS -l walltime=2:00:00
#
cd ${HOME}/Delphes-3.0.12/
if [ -f ${rootname} ]
then
	rm ${rootname}
fi
${HOME}/Delphes-3.0.12/DelphesHepMC ${HOME}/Delphes-3.0.12/examples/delphes_card_CMS_cleanTk.tcl ${rootname} ${hepmcname}