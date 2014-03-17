#!/bin/bash

# Submits lots of qcd making jobs

for i in {21..40}
do
	echo "Doing qcdb_$i.hepmc"
	qsub -v name="qcdb_$i.hepmc" runQCDbBatch.sh
	sleep 10
done
