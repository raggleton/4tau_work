#!/bin/bash

for f in qcdJob.sh.o*
do
	process=$(cat $f |head -n 4 | tail -n 1 | awk '{print $2;}')
	echo $process
	number=${f#qcdJob.sh.o}
	echo $number
	if [ "$process" == "qcdb" ]
	then
		echo "qcdb"
		mv $f ../QCDb_HLT/jobOutput/
		mv qcdJob.sh.e$number ../QCDb_HLT/jobOutput/
	elif [ "$process" == "qcdscatter" ]
	then	
		echo "qcdscatter"
		mv $f ../QCDScatter_HLT/jobOutput/
		mv qcdJob.sh.e$number ../QCDScatter_HLT/jobOutput/
	fi	
done
