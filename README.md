4tau_work
=========

For Robin's work on NMSSM H1 -> 2a1 -> 4tau

Monte Carlo stuff is in `pythia8` and `Delphes3` folders.

Installing & compiling against Boost
====================================

http://www.boost.org/doc/libs/1_55_0/more/getting_started/unix-variants.html

1. Download `boost_1_55_0.tar.bz2`.
2. Extract `tar --bzip2 -xf /path/to/boost_1_55_0.tar.bz2` (say to `~/boost_1_55_0`).
3. Now need to build the separate Boost libraries 
	>The only Boost libraries that must be built separately are:
	>
	>Boost.Chrono
	>Boost.Context
	>Boost.Filesystem
	>Boost.GraphParallel
	>Boost.IOStreams
	>Boost.Locale
	>Boost.MPI
	>Boost.ProgramOptions
	>Boost.Python (see the Boost.Python build documentation before building and installing it)
	>Boost.Regex
	>Boost.Serialization
	>Boost.Signals
	>Boost.System
	>Boost.Thread
	>Boost.Timer
	>Boost.Wave

4. Create installation folder for those libraries: `mkdir ~/boost_1_55_0_install`
5. In the original boost folder, do `./bootsrap.sh --prefix=$(HOME)/boost_1_55_0_install`
6. Then run `./b2 install`. This should install libraries to  `~/boost_1_55_0_install/lib`
7. In the Delphes Makefile, change the following:
	- Add `-I $(HOME)/boost_1_55_0 -I $(HOME)/boost_1_55_0_install/include` to `CCXFLAGS` (line 18)
	- Add `-L/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib -lboost_program_options` to `DELPHES_LIBS` (line 19)
	- (May also require adding the lib folder to LD_LIBRARY_PATH, but even if you do, you still need the -L with path so ???)
8. Note that if you use another boost library (like Filesystem) you'll need to add `-lboost_XXX` to `DELPHES_LIBS`
9. For Boost Program Options, you can now follow the tutorial. 