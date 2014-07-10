4tau_work
=========

For Robin's work on NMSSM H1 -> 2a1 -> 4tau

Monte Carlo stuff is in `pythia8` and `Delphes3` folders.
You'll need Boost both Pythia and Delphes.

Sorry in advance if this is a bit rubbish - it's constantly in a state of flux, and no-one but me uses it. (Standard weak HEP excuse)

## Setup BOOST 

- See below

## Setup PYTHIA8

- Download tarball, extract
- Follow README.HepMC

## Setup Delphes3

** Requires version 3.1.2 or greater** (basically need track class to have impact parameter variables Dxy,Xd, Yd, etc - check in `DelphesClasses.h`)

- Download tarball, extract
- `cd` into folder, run `make`
- Edit and run `setupDelphesLinks.sh` to create symlinks to files in your Delphes folder

## Workflow

**Signal**:
- Make signal MC using CalcHEP. Outputs as LHE file.
- Hadronise it using Pythia: 
```
make mainLHEhadronise
./mainLHEhadronise.exe # to show usage
./mainLHEhadronise.exe <LHEfilename & path> <output HEPMC filename>
```
This outputs a HEPMC file.

**QCD[b/c/q-g scatter]**:
- Make it all in Pythia using one executable
```
make mainQCD
./mainQCD.exe # to show usage
```
Once you have a HEPMC file, you now run it through Delphes to produce a ROOT file:
- Run `./DelphesHepMC` for relevant options.
- Generally, you'll want to use `DelphesCards\delphes_card_cms_bare.tcl`. If you're using a version of Delphes older than 3.1.2, use `delphes_card_CMS_bare_OLD.tcl`.
- Use `Scripts/submitDelphesLots.sh` or `Scripts/submitDelphesSingle.sh` for use on PBS batch system, or `Scripts/runDelphesLots.sh` for lots of local jobs.

## Delphes Analysis Programs:
- `basicScript`: template script for making a new program (see below) **OUT OF DATE**
- `mainAnalysis`: does lots of plots, like pT, track eta Vs phi, soft track distributions. Used for QCDb rejection studies **OUT OF DATE**
- `massPlots`: makes plots of invariant masses of mu+tk, and calculates & plots correlation coefficients for backgrounds
- `IP`: plots impact parameters of things. Not really kept up to date. **OUT OF DATE**

## Making a new Delphes analysis program: **OUT OF DATE**
- Look at `Scripts/createNewScript.sh` - you'll need to edit the relevant paths
- Run `createNewScript.sh <myscriptName>`. Make sure it puts the new scripts in your `Delphes/examples` folder
- Go to your Delphes directory, run `./configure` to pick up new script
- Look at `Scripts/createMakefile.sh`, modify paths as necessary
- Run `Scripts/createMakefile.sh` to modify the Delphes makefile to add support for C++11 and Boost.
- Now run `make`
- To run your program, do `./runMyCoolNewProgram` in the main Delphes installation folder.

You can now edit your new analysis code. To re-make it, just run `makeScripts`. This will re-make ALL the analysis programs. See `makeMain` and `makeMass` for examples on how to write a script to re-make just one program.

### Notes

- There are old Pythia programs, `mainQCDb`, `mainQCDScatter`, etc. **These are are all old, and buggy. Do not use them.** Only use `mainQCD` (plus, it does everything with command-line options, so win).


Installing & compiling against Boost
====================================

http://www.boost.org/doc/libs/1_55_0/more/getting_started/unix-variants.html

1. Download `boost_1_55_0.tar.bz2`.
2. Extract `tar --bzip2 -xf /path/to/boost_1_55_0.tar.bz2` (say to `~/boost_1_55_0`).
3. Now need to build the separate Boost libraries 
	>- The only Boost libraries that must be built separately are:
	>
	>
	>- Boost.Chrono
	>- Boost.Context
	>- Boost.Filesystem
	>- Boost.GraphParallel
	>- Boost.IOStreams
	>- Boost.Locale
	>- Boost.MPI
	>- Boost.ProgramOptions
	>- Boost.Python (see the Boost.Python build documentation before building and installing it)
	>- Boost.Regex
	>- Boost.Serialization
	>- Boost.Signals
	>- Boost.System
	>- Boost.Thread
	>- Boost.Timer
	>- Boost.Wave

4. Create installation folder for those libraries: `mkdir ~/boost_1_55_0_install`
5. In the original boost folder, do `./bootsrap.sh --prefix=$(HOME)/boost_1_55_0_install`
6. Then run `./b2 install`. This should install libraries to  `~/boost_1_55_0_install/lib`
7. In the Delphes Makefile, change the following:
	- Add `-I $(HOME)/boost_1_55_0 -I $(HOME)/boost_1_55_0_install/include` to `CCXFLAGS` (line 18)
	- Add `-L/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib -lboost_program_options` to `DELPHES_LIBS` (line 19)
	- Also adding the `lib` folder to `LD_LIBRARY_PATH` in your `~/.bashrc`:
		```
		export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib/
		```
8. Note that if you use another boost library (like Filesystem) you'll need to add `-lboost_XXX` to `DELPHES_LIBS`
9. For Boost Program Options, you can now follow the tutorial. 

```
LD_LIBRARY_PATH=/cm/shared/apps/gcc/4.7.0/lib:/cm/shared/apps/gcc/4.7.0/lib64:/cm/shared/languages/Python-2.7.6/lib:/cm/shared/apps/torque/4.2.4.1/lib:/cm/shared/apps/moab/7.2.2/lib:/cm/shared/tools/git-1.8.4.2/lib:/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib/:/panfs/panasas01/phys/ra12451/root/root/lib/
```

## Analysis Note

First time setup:
```
svn co -N svn+ssh://svn.cern.ch/reps/tdr2 myDir
# svn co -N svn+ssh://raggleto@svn.cern.ch/reps/tdr2 myDir
cd myDir
svn update utils
svn update -N notes
svn update notes/AN-14-007
eval `notes/tdr runtime -sh` # (for bash. use -csh for tcsh.)
cd notes/AN-14-007/trunk/
tdr --style=an b AN-14-007
```

To update/remake:
```
cd myDir
svn update notes/AN-14-007
eval `notes/tdr runtime -sh` # (for bash. use -csh for tcsh.)
cd notes/AN-14-007/trunk/
tdr --style=an b AN-14-007
```

See https://twiki.cern.ch/twiki/bin/view/CMS/Internal/TdrProcessing