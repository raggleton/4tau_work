#!/bin/bash

for i in {1..10}
do
	cp nmssm_4tau_ma8.txt nmssm_4tau_ma8_$i.txt
	sed -i "s/NUMEVT/$i/g" nmssm_4tau_ma8_$i.txt   
	echo "Doing nmssm_4tau_ma8_${i}.txt"
    if [ $i -eq 1 ]
    then
        submitString=`qsub -N nmssm_${i} -v filename="nmssm_4tau_ma8_${i}.txt" runCalchepBatch.sh`
        echo "qsub -N nmssm_${i} -v filename="nmssm_4tau_ma8_${i}.txt" runCalchepBatch.sh"
        echo "$submitString"
        jobID=`echo $submitString | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
    else
        submitString=`qsub  -W depend=afterok:${jobID} -N nmssm_${i} -v filename="nmssm_4tau_ma8_${i}.txt" runCalchepBatch.sh`
        echo "qsub  -W depend=afterok:${jobID} -N nmssm_${i} -v filename="nmssm_4tau_ma8_${i}.txt" runCalchepBatch.sh"
        echo "$submitString"
    	jobID=`echo $submitString | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
    fi
done
