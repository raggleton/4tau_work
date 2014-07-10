#!/bin/bash


# Edit these as necessary
SOURCE="/panfs/panasas01/phys/ra12451/4Tau/4tau_work/Delphes3/"
DEST_EXAMPLES="/panfs/panasas01/phys/ra12451/Delphes-3.1.2/examples/"
DEST_MAIN="/panfs/panasas01/phys/ra12451/Delphes-3.1.2/"

cd -- "$SOURCE"
find -maxdepth 1 -type f -exec ln -s "$SOURCE"/{} "$DEST_EXAMPLES"/{} \;
# find . -maxdepth 1 -name make* -type f -exec ln -s "$SOURCE"/{} "$DEST_MAIN"/{} \;

SOURCE2="$SOURCE/Scripts"
echo $SOURCE2
cd -- "$SOURCE2"
find -maxdepth 1 -type f -exec ln -s "$SOURCE2"/{} "$DEST_EXAMPLES"/{} \;
# find "$SOURCE2" -maxdepth 1 -name make* -type f -exec ln -s "$SOURCE2"/{} "$DEST_MAIN"/{} \;

SOURCE2="$SOURCE/DelphesCards"
cd -- "$SOURCE2"
find -maxdepth 1 -type f -exec ln -s "$SOURCE2"/{} "$DEST_EXAMPLES"/{} \;

# A few custom ones
mv "$DEST_EXAMPLES/createMakefile.sh" "$DEST_MAIN"
mv "$DEST_EXAMPLES/"make* "$DEST_MAIN"
# mv "$DEST_EXAMPLES/Makefile" "$DEST_MAIN"
rm "$DEST_EXAMPLES/Makefile"
exit