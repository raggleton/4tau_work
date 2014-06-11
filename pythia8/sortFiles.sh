#!/bin/bash

for f in $(ls /panfs/panasas01/phys/ra12451/pythia8183/examples/jobOutput/*)
do
	# echo $f
	mline=$(head -n 2 $f | tail -n 1)
	if echo $mline| grep "qcdScatter"
	then
		echo $mline
		echo $f
		mv $f /panfs/panasas01/phys/ra12451/pythia8183/examples/jobOutput/QCDScatter/
	fi

	# if echo $mline | grep "qcdb_pthatmin20_Mu17_Mu8_"
	# 	then
	# 	echo $mline
	# 	echo $f
	# fi
done