#!/bin/bash

# Script that converts the event listing from Pythia 
# into a Graphviz file to be plotted with dot
# e.g. 
# ./makeGraphFromPythia > myEvent.gv
# dot -Tpng myEvent.gv -o myEvent.png

# Filename for output file
outputFile="myEvent.gv"

# array to hold names of particle for given No in event listing
nameArr=($(awk '{print $1":"$3}' testLine.txt))
# Store index of first daughter
d1Arr=($(awk '{print $7}' testLine.txt))
# Store index of second daughter
d2Arr=($(awk '{print $8}' testLine.txt))

nParticles=${#nameArr[@]}
let "nParticles--"
echo "digraph g {" > $outputFile
echo "    rankdir = LR;" >> $outputFile
for i in $(eval echo {0..$nParticles})
# for i in {2..10}
do
	# Create string to be written to fle
	gEntry="    \"${nameArr[$i]}\" -> "
	
	finalState=false
	
	# Check if d2==0 (we dont want all particles!)
	if [[ ${d2Arr[$i]} -eq 0 ]]
	then
		# Check if both daughters are == 0 - then final state
		if [[ ${d1Arr[$i]} -eq 0 ]]
		then
			finalState=true
		else
			d1No=${d1Arr[$i]}
			dau=${nameArr[$d1No]}
			gEntry+='"$dau"'
		fi
	else
		gEntry+="{ "
		d1=${d1Arr[$i]}
		d2=${d2Arr[$i]}
		for d in $(eval echo {$d1..$d2})
		do
			dau=${nameArr[$d]}
			gEntry+="\"$dau\" "
		done
		gEntry+=" }"
	fi
	
	if [[ !finalState ]]; then
		echo "$gEntry"  >> $outputFile
	fi
done
echo "}"  >> $outputFile
