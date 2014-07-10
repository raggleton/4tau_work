#!/bin/bash


# Edit these as necessary
SOURCE="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/"
DEST="/panfs/panasas01/phys/ra12451/Delphes-3.1.2/examples/"

cd -- "$SOURCE"
find -maxdepth 1 -type f -exec ln -s "$SOURCE"/{} "$DEST"/{} \;

SOURCE="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/DelphesCards"
cd -- "$SOURCE"
find -maxdepth 1 -type f -exec ln -s "$SOURCE"/{} "$DEST"/{} \;

SOURCE="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/Scripts"
cd -- "$SOURCE"
find -maxdepth 1 -type f -exec ln -s "$SOURCE"/{} "$DEST"/{} \;

exit