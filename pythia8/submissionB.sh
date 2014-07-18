#!/bin/bash

#for i in {1..10}
#do
#  	qsub -v args="-n 5000 -p qcdb --writeHLT --name b_5000_${i} --notMuOnly" qcdJob.sh
#	echo "qsub -v args="-n 5000 -p qcdb --writeHLT --name b_5000_${i} --notMuOnly" qcdJob.sh"
#        sleep 10
#done

for i in {602..800}
do
        echo "qsub -N b$i -v args=\"-n 5000 -p qcdb --writeHLT --name b_5000_${i} --seed ${i}\" qcdJob.sh"
	qsub -N b$i -v args="-n 5000 -p qcdb --writeHLT --name b_5000_${i} --seed ${i}" qcdJob.sh
	sleep 10
done

