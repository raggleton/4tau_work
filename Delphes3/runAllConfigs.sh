#!/bin/bash

# Run ALL configurations: (sig,QCDb)*(muRand, mu Ordered)*(HLT, NoHLT)
# FILE="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/mainAnalysis.C"
FILE="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/test.txt"

printConfig(){
	head -22 $FILE | tail -4
}

makeAnalysis(){
	touch examples/mainAnalysis.cpp && make
}

checkPrevJob(){
	# checks to see if previous job finished
	# -q to supress output
	while qstat -u ra12451 | grep runAnalysis.sh
	do
		date +"%T"
		echo "Job still running..."
		sleep 30m
	done
}

submitJob(){
	qsub runAnalysis.sh
}

printTime(){
	date +"%T"
}

# Make sure nothing is running first!
checkPrevJob

########
# SIGNAL
########
sed -i.bak 's:doSignal = .*; :doSignal = true\; :' test.txt
	
	######
	# HLT
	######
	sed -i.bak 's:doHLT = .*; :doHLT = true\; :' test.txt

		#########
		# muRand
		#########
		date +"%T"
		echo "SIGNAL-HLT-MU_RAND"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = true\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob

		############
		# muOrdered 
		############
		date +"%T"
		echo "SIGNAL-HLT-MU_ORDERED"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = false\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob

	########
	# NoHLT
	########
	sed -i.bak 's:doHLT = .*; :doHLT = false\; :' test.txt
	
		#########
		# muRand
		#########
		date +"%T"
		echo "SIGNAL-NOHLT-MU_RAND"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = true\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob

		############
		# muOrdered 
		############
		date +"%T"
		echo "SIGNAL-NOHLT-MU_ORDERED"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = false\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob

#########
# QCDb_mu
########
sed -i.bak 's:doSignal = .*; :doSignal = false\; :' test.txt
	
	######
	# HLT
	######
	sed -i.bak 's:doHLT = .*; :doHLT = true\; :' test.txt

		#########
		# muRand
		#########
		date +"%T"
		echo "QCD-HLT-MU_RAND"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = true\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob

		############
		# muOrdered 
		############
		date +"%T"
		echo "QCD-HLT-MU_ORDERED"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = false\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob

	########
	# NoHLT
	########
	sed -i.bak 's:doHLT = .*; :doHLT = false\; :' test.txt
	
		#########
		# muRand
		#########
		date +"%T"
		echo "QCD-NOHLT-MU_RAND"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = true\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob

		############
		# muOrdered 
		############
		date +"%T"
		echo "QCD-NOHLT-MU_ORDERED"
		sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = false\; :' test.txt
		printConfig
		makeAnalysis
		submitJob
		checkPrevJob
