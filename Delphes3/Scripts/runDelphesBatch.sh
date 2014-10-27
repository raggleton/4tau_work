#!/bin/bash
#
#PBS -l walltime=3:00:00
#
cd ${HOME}/Delphes-3.1.2/
if [ -f ${rootname} ]
then
	rm ${rootname}
fi
${HOME}/Delphes-3.1.2/DelphesHepMC ${HOME}/Delphes-3.1.2/examples/delphes_card_CMS_bare.tcl ${rootname} ${hepmcname}
