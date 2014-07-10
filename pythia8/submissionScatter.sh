#!/bin/bash

for i in {1..200}
do
	echo "qsub -N scatter$i -v args="-n 250 -p qcdscatter --writeHLT --name scatterbc_250_${i}" qcdJob.sh"
  	qsub -N scatter$i -v args="-n 250 -p qcdscatter --writeHLT --name scatterbc_250_${i}" qcdJob.sh
        sleep 10
done
