#!/bin/bash

for name in ${HOME}/pythia8183/Signal_1prong_500K_bare/GG_H_aa_8_4taus_decay_500K_*.hepmc
do
	# Pure filename, no path
	filename=${name##*/}
	# echo $filename

	# Get whether NoHLT or HLT
	hltConf=""
	if [[ $filename == *NoHLT* ]]
		then
		hltConf="NoHLT"
	else
		hltConf="HLT"
	fi

	# Get the file number
	tmp=${filename#GG_H_aa_8_4taus_decay_500K_}
	num=${tmp%-single*}
	# echo $num

	# Set filename
	# Set which Delphes config used
	delphConf="bare"
	echo "Doing "$filename
	out="signal_1prong_500K_${num}_${hltConf}_${delphConf}"	
	echo $out
	qsub -v hepmcname=${name},rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/Signal_1prong_500K_bare/"${out}.root runDelphesBatch.sh
	sleep 10
done
