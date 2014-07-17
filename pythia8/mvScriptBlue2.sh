#!/bin/bash

# Move output of jobs into folders
# Assumes that job output logs from PBS are named something like b789.o123456
# where the "789" bit corresponds to the hepMC file number and random generator seed
filestem="b[0-9]*.o"
for f in $filestem*
do
	process=$(cat $f |head -n 4 | tail -n 1 | awk '{print $2;}')
	# echo $process
	
	# Get PBS job ID number from line 4
	number=${f#$filestem}
	# echo $number
	
	# Get the file number (the bit before the .oxxxxx)
	fnumber=${f%.o*}
	fnumber=${fnumber#b}
	# echo $fnumber

	# Get the hepmc filenames from line 11
	file1=$(cat $f |head -n 11 | tail -n 1 | awk '{print $3;}')
	file2=$(cat $f |head -n 11 | tail -n 1 | awk '{print $5;}')
	# alternatively could use fnumber from above...

	folder=""
	if [ "$process" == "qcdb" ]
	then
		# echo "qcdb"
		folder="QCDb_HLT"
	elif [ "$process" == "qcdscatter" ]
	then	
		# echo "qcdscatter"
		folder="QCDScatter_HLT"
	fi	
	
	# Check folder actually exists first
	if [ ! -d "../$folder" ]; then
		mkdir ../"$folder"
		mkdir ../"$folder"/jobOutput
	fi

	echo "mv $f ../"$folder"/jobOutput/$f"
	echo "mv ${f/o/e} ../"$folder"/jobOutput/${f/o/e}"
	echo "mv b_5000_${fnumber}_progress.txt ../"$folder"/jobOutput/b_5000_${fnumber}_progress.txt"
	echo "mv $file1 ../"$folder"/$file1"
	echo "mv $file2 ../"$folder"/$file2"
done
