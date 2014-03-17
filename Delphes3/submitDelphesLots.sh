#!/bin/bash

for name in ${HOME}/pythia8183/examples/QCDb_*.hepmc
do
	filename=${name#*examples/}
	echo "Doing "$filename
	num=${filename#QCDb_}
	out="QCDb_mu_$num"	
	qsub -v hepmcname=${name},rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/examples/QCDb_mu_cleanTk/"${out%hepmc}root runDelphesBatch.sh
	sleep 10
done