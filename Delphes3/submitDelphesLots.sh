#!/bin/bash

# for name in ${HOME}/pythia8183/examples/qcdb_*.hepmc
# for name in ${HOME}/pythia8183/QCDb_mu_pthatmin20/qcdb_pthatmin20_*.hepmc
for name in ${HOME}/pythia8183/QCDb_mu_pthatmin20_Mu17_Mu8/qcdb_pthatmin20_Mu17_Mu8_*.hepmc
do
	filename=${name##*/}
	num=${filename##*_}
	pureNum=${num%.hepmc}
	# echo $pureNum
	if [ "$pureNum" -gt "120" ]
	then
		echo "Doing "$filename
		out="QCDb_mu_pthatmin20_Mu17_Mu8_$pureNum"	
		# echo $out
		# echo $name
		qsub -v hepmcname=${name},rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDb_mu_pthatmin20_Mu17_Mu8_bare/"${out}.root runDelphesBatch.sh
		sleep 10
	fi
done
