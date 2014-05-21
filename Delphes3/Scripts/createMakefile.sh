#!/bin/bash

# This modifies the Delphes makefile to add in BOOST, ROOT stuff, C++11 support
# Run this any time you run ./configure as it'll overwrite it
# 
CXXMOD="-I \$(HOME)/boost_1_55_0 -I \$(HOME)/boost_1_55_0_install/include -std=c++11"
LIBMOD="-L/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib -lboost_program_options -lboost_filesystem"

echo "Appending to Makefile..."
echo "To CXXFLAGS: " $CXXMOD
echo "To DELPHES_LIBS: " $LIBMOD

# This gets the first CXXFLAGS and DELPHES_LIBS lines 9had to manually restrict line number range todo this)
# The & is the matched string (so CXXFLAGS += ...)
# The '"$MY_VAR"' is necessary to ensure spaces etc are dealt with properly
# sed is bloody horrible!
sed -i -e '17,21 s:^CXXFLAGS .*:& '"$CXXMOD"':' -e '17,21 s:^DELPHES_LIBS .*:& '"$LIBMOD"':' Makefile
