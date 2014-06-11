#!/usr/bin/python

# To add total number generated
# Get from the output files
# import os, os.path
import glob
import re

def grep(filename, searchString):
	returnFilename=""
	with open(filename,"r") as file:
		for line in file:
			if re.search(re.escape(searchString), line):
				returnFilename = filename
	return returnFilename


# Do QCD scatter first
# totalScatter = 0.0
# listing = glob.glob("/panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDScatter_mu_pthatmin20_Mu17_Mu8_bare/QCDScatter_mu_pthatmin20_Mu17_Mu8_*.root")
# nfilesScatter = len(listing)
# scatterfilesCounted = 0.0 # count number that actually contribute to sum
# print nfilesScatter, " QCD Scatter files"

# # # Loop over all ROOT files
# # for i in range (1,nfilesScatter+1):
# outputfiles = glob.glob("/panfs/panasas01/phys/ra12451/pythia8183/examples/jobOutput/QCDScatter/*")
# # for i in range (1,3):
# for i in range (1,nfilesScatter+1):
# 	# print "QCDScatter_mu_pthatmin20_Mu17_Mu8_"+str(i)+": "
# 	# For each, loop over all jobOutput files, see if any contain the hepmc filename that corresponds to the ROOT filename
# 	numberEvents = 0
# 	print len(outputfiles),"files to look over"
	
# 	for f in outputfiles:
# 		searchterm = "qcdScatter_pthatmin20_"+str(i)+"_HLT"
# 		fileResult = grep(f,searchterm)
# 		if fileResult:
# 			with open(fileResult,"r") as file:
# 				for line in file:
# 					if re.search("\| sum",line):
# 						# print line
# 						values = line.split()
# 						numberEvents = int(values[5])
# 						totalScatter += int(values[5])
# 						if int(values[5]) > 0:
# 							scatterfilesCounted += 1
# 						break

# 			break

# 	if fileResult:
# 		outputfiles.remove(fileResult)
# 	print "QCDScatter_mu_pthatmin20_Mu17_Mu8_"+str(i)+": "+str(numberEvents)
# 	print fileResult

# print "TOTAL NUMBER OF SCATTER EVENTS: ", totalScatter
# print "SHOULD BE: ", nfilesScatter*(totalScatter/float(scatterfilesCounted))
# print "AVERAGE EVT/FILE: ", (totalScatter/float(nfilesScatter))
# print "FILES COUNTED: ", scatterfilesCounted

# Do QCDb next
totalB = 0.0
listing = glob.glob("/panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDb_mu_pthatmin20_Mu17_Mu8_bare/QCDb_mu_pthatmin20_Mu17_Mu8_*.root")
nfilesB = len(listing)
BfilesCounted = 0.0 # count number that actually contribute to sum
print nfilesB,"QCDb files"

# Loop over all ROOT files
outputfiles = glob.glob("/panfs/panasas01/phys/ra12451/pythia8183/examples/jobOutput/runQCDb*")
for i in range (1,nfilesB+1):
	print "QCDb_mu_pthatmin20_Mu17_Mu8_"+str(i)+": "
	# For each, loop over all jobOutput files, see if any contain the hepmc filename that corresponds to the ROOT filename
	numberEvents = 0
	print len(outputfiles),"files to look over"
	
	for f in outputfiles:
		search = "qcdb_pthatmin20_"+str(i)
		fileResult = grep(f,search)
		if fileResult:
			print fileResult
			with open(fileResult,"r") as file:
				for line in file:
					if re.search("\| sum",line):
						# print line
						values = line.split()
						numberEvents = int(values[5])
						totalB += int(values[5])
						if int(values[5]) > 0:
							BfilesCounted += 1
						break
			break
	if fileResult:
		outputfiles.remove(fileResult)
	print "QCDb_mu_pthatmin20_Mu17_Mu8_"+str(i)+": "+str(numberEvents)
	print fileResult

print "TOTAL NUMBER OF B EVENTS COUNTED: ", totalB
print "SHOULD BE: ", nfilesB*(totalB/float(BfilesCounted))
print "AVERAGE EVT/FILE: ", (totalB/float(BfilesCounted))
print "FILES COUNTED: ", BfilesCounted
