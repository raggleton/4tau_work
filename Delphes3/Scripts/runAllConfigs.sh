#!/bin/bash

#######################################################################
# Run ALL configurations: (sig,QCDb)*(muRand, mu Ordered)*(HLT, NoHLT)
#######################################################################

FILE="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/mainAnalysis.C"
# FILE="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/$FILE"

printConfig(){
	# Print the bool switches to check
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

printMakeSubmitCheck(){
	printConfig
	makeAnalysis
	submitJob
	checkPrevJob
}

# Set switches in mainAnalysis.C
doSignal(){
	sed -i.bak 's:doSignal = .*; :doSignal = true\; :' $FILE
}

doQCD(){
	sed -i.bak 's:doSignal = .*; :doSignal = false\; :' $FILE
}

doHLT(){
	sed -i.bak 's:doHLT = .*; :doHLT = true\; :' $FILE
}

doNoHLT(){
	sed -i.bak 's:doHLT = .*; :doHLT = false\; :' $FILE
}

doMuRand(){
	sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = true\; :' $FILE
}

doMuOrdered(){
	sed -i.bak 's:swapMuRandomly = .*; :swapMuRandomly = false\; :' $FILE
}

###########################
# MAIN PROGRAM STARTS HERE
###########################

# Make sure nothing is running first!
checkPrevJob

########
# SIGNAL
########
doSignal	
	######
	# HLT
	######
	doHLT
		#########
		# muRand
		#########
		printTime
		echo "SIGNAL-HLT-MU_RAND"
		doMuRand
		printMakeSubmitCheck

		############
		# muOrdered 
		############
		printTime
		echo "SIGNAL-HLT-MU_ORDERED"
		doMuOrdered
		printMakeSubmitCheck

	########
	# NoHLT
	########
	doNoHLT
	
		#########
		# muRand
		#########
		printTime
		echo "SIGNAL-NOHLT-MU_RAND"
		doMuRand
		printMakeSubmitCheck

		############
		# muOrdered 
		############
		printTime
		echo "SIGNAL-NOHLT-MU_ORDERED"
		doMuOrdered
		printMakeSubmitCheck

#########
# QCDb_mu
#########
doQCD

	######
	# HLT
	######
	doHLT

		#########
		# muRand
		#########
		printTime
		echo "QCD-HLT-MU_RAND"
		doMuRand
		printMakeSubmitCheck

		############
		# muOrdered 
		############
		printTime
		echo "QCD-HLT-MU_ORDERED"
		doMuOrdered
		printMakeSubmitCheck

	########
	# NoHLT
	########
	doNoHLT
	
		#########
		# muRand
		#########
		printTime
		echo "QCD-NOHLT-MU_RAND"
		doMuRand
		printMakeSubmitCheck

		############
		# muOrdered 
		############
		printTime
		echo "QCD-NOHLT-MU_ORDERED"
		doMuOrdered
		printMakeSubmitCheck
