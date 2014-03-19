#!/bin/bash

# This script makes a new Delphes analysis script based off basicScript.C
#
# Argument: new script name

if [ $# -ne 1 ]
then
	echo "Need to specify new script name!"
	exit
fi

top="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/"

cp $top/basicScript.C $top$1.C
cp $top/basicScript.cpp $top$1.cpp

ln -s $top$1.C /panfs/panasas01/phys/ra12451/Delphes-3.0.12/examples/$1.C
ln -s $top$1.cpp /panfs/panasas01/phys/ra12451/Delphes-3.0.12/examples/$1.cpp

sed -i "s/basicScript/$1/g" $top$1.C
sed -i "s/basicScript/$1/g" $top$1.cpp
