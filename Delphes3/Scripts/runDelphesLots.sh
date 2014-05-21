#!/bin/bash

# Script to do lots of Delphes submissions on pbs batch system

# for name in ${HOME}/pythia8183/examples/qcdb_*.hepmc
# for name in ${HOME}/pythia8183/QCDb_mu_pthatmin20/qcdb_pthatmin20_*.hepmc
# for name in ${HOME}/pythia8183/QCDb_mu_pthatmin20_Mu17_Mu8/qcdb_pthatmin20_Mu17_Mu8_*.hepmc
for name in ${HOME}/pythia8183/QCDc_mu_pthatmin20_Mu17_Mu8/qcdc_pthatmin20_*.hepmc
# for name in ${HOME}/pythia8183/QCDb_mu_pthatmin20/qcdb_pthatmin20_*.hepmc
do
	filename=${name##*/}
	# num=${filename##*_}
	# pureNum=${num%.hepmc}
	# stem=${filename%_NoHLT.hepmc}
	stem=${filename%_HLT.hepmc}
	# pureNum=${stem#qcdb_pthatmin20_Mu17_Mu8_}
	pureNum=${stem#qcdc_pthatmin20_}
	if [ "$pureNum" -gt "100" ]
	then
		echo "Doing "$filename
		echo "Writing to "${out}.root
		out="QCDc_mu_pthatmin20_Mu17_Mu8_$pureNum"	
		# echo $out
		# echo $name
		# qsub -v hepmcname=${name},rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDb_mu_pthatmin20_Mu17_Mu8_bare/"${out}.root runDelphesBatch.sh
		# qsub -v hepmcname=${name},rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDc_mu_pthatmin20_Mu17_Mu8_bare/"${out}.root runDelphesBatch.sh
		# sleep 10
		../DelphesHepMC ../examples/delphes_card_CMS_bare.tcl /panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDc_mu_pthatmin20_Mu17_Mu8_bare/${out}.root $name
	fi
done
