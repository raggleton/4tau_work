#!/bin/bash

# This script makes a new Delphes analysis script based off basicScript.C.
# Run from your github dir, so the new scripts are kept in git dir, 
# and symlinks put in DELPHESDIR. Then you can commit changes.

# BEFORE RUNNING: change DELPHESDIR
# Argument: new script name

if [ $# -ne 1 ]
then
	echo "Need to specify new script name!"
	exit
fi

DELPHESDIR=${HOME}/Delphes-3.1.2
# top="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/"
GITDELPHESDIR=$(pwd)
GITDELPHESDIR=${GITDELPHESDIR%Scripts}
echo DELPHESDIR=$DELPHESDIR
echo GITDELPHESDIR=$GITDELPHESDIR
cp -i $GITDELPHESDIR/basicScript.C $GITDELPHESDIR/$1.C
cp -i $GITDELPHESDIR/basicScript.cpp $GITDELPHESDIR/$1.cpp

ln -sf $GITDELPHESDIR/$1.C $DELPHESDIR/examples/$1.C
ln -sf $GITDELPHESDIR/$1.cpp $DELPHESDIR/examples/$1.cpp


if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform        
    # sed on osx requires a backup extension. here I've told it to ignore it
    sed -i '' s/basicScript/$1/g $GITDELPHESDIR/$1.C
    sed -i '' s/basicScript/$1/g $GITDELPHESDIR/$1.cpp
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Do something under Linux platform
    sed -i s/basicScript/$1/g $GITDELPHESDIR/$1.C
    sed -i s/basicScript/$1/g $GITDELPHESDIR/$1.cpp
fi
