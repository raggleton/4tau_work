#ifndef PROGRAMOPTS_H
#define PROGRAMOPTS_H

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <exception>

// ROOT headers
// #include "TStyle.h"

// BOOST headers
// Need to add
// -I $(HOME)/boost_1_55_0 -I $(HOME)/boost_1_55_0_install/include
// to CXXFLAGS in Delphes/Makefile
// and
// -L/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib -lboost_program_options -lboost_filesystem
// to DELPHES_LIBS in Delphes/Makefile
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using std::cout;
using std::endl;

namespace fs = boost::filesystem;
namespace po = boost::program_options;

/**
 * This header contains all the code relating to the program options
 * for Delphes analysis scripts (ie adding of ROOT files for analysis)
 *
 * Robin Aggleton 2014
 */

// for rescaling deltaR(mu-tk)
// get this value from running reweightingVariables program, then running muonReweighting.C script on output of that
// then put in the two numbers here
// alpha = <deltaR in data> / <deltaR in MC>
// where < > indicate mean. So alpha should be < 1
const double alpha = 0.150897/0.175266;

// Global enum for the MC source
// Note convention of all lower case
enum MCsource { signal, qcdb, qcdc, qcdscatter, qcdall, test };

// Have to define << and >> ops to get enum to work with boost::lexical_cast
// and program_options, and also so we can cout the enum easily
std::ostream& operator<<(std::ostream& out, const MCsource value) {

    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value) {
        PROCESS_VAL(signal);
        PROCESS_VAL(qcdb);
        PROCESS_VAL(qcdc);
        PROCESS_VAL(qcdscatter);
        PROCESS_VAL(qcdall);
        PROCESS_VAL(test);
    }
#undef PROCESS_VAL

    return out << s;
}

// Must be a way to improve this using static_cast
std::istream & operator>>(std::istream & in, MCsource & value) {
  std::string token;
  if (in >> token) {
  	if (token == "signal")
  		value = signal;
  	else if (token == "qcdb")
  		value = qcdb;
  	else if (token == "qcdc")
  		value = qcdc;
  	else if (token == "qcdscatter")
  		value = qcdscatter;
  	else if (token == "qcdall")
  		value = qcdall;
  	else if (token == "test")
  		value = test;
  	else
  		throw runtime_error("Invalid string cast to enum");
  }
  return in;
}

/**
 * This class is to handle program options using boost::program_options
 */
class ProgramOpts
{
	private:
		MCsource source; // do signal or qcd(b)(c)(scatter) or test (a qcdb file)
		bool doSignal; // do signal, or not signal (some parts are signal only)
		bool doMu; // for QCDb - either inclusive decays or mu only decays - DEPRECIATED
		bool swapMuRandomly; // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
		bool doHLT; // whether to use MC that has HLT cuts already applied.
		int  nEvents; // how many events to run over, -1 = all events
		bool verbose; // output debugging info to screen
		double dR; // deltaR(mu-mu) cut
		bool doRescale; // whether to rescale the track eta & phi to match with data
		double rescale; // rescale value for rescaling the track eta & phi. 1 = no rescaling

	public:
		// constructor, parses input
		ProgramOpts(int argc, char* argv[]):
			// some sensible defaults
			source(signal),
			doSignal(true),
			doMu(true),
			swapMuRandomly(true),
			doHLT(true),
			nEvents(-1),
			verbose(false),
			dR(2.),
			doRescale(false),
			rescale(1.0)
		{
			po::options_description desc("\nAllowed options");
			desc.add_options()
				("help,h", "Produce help message")
				("source,s", po::value<MCsource>(&source), "Process to run: signal [default], qcdb, qcdc, qcdscatter, test (a qcdb file)")
				("swapMuRandomly", po::value<bool>(&swapMuRandomly), "TRUE [default] - mu 1,2 randomly assigned, FALSE - mu 1,2 pT ordered")
				("doHLT", po::value<bool>(&doHLT), "TRUE [default] - use samples with HLT_Mu17_Mu8 during generation, FALSE - no HLT cuts")
				("number,n", po::value<int>(&nEvents), "Number of events to run over. -1 for all [default]")
				("dR", po::value<double>(&dR), "Set deltaR(mu-u) cut value (default = 2).")
				("verbose,v", "Output debugging statements")
				("doRescale", "Whether to rescale track eta & phi to match data. FALSE by default.")
			;

			po::variables_map vm;
			try {
				po::store(po::parse_command_line(argc, argv, desc), vm);
			} catch (boost::program_options::invalid_option_value e) {
				cout << "Invalid option value: " << e.what() << endl;
				cout << desc << endl;
				cout << "Exiting" << endl;
				exit(1); // NOT ELEGANT - DO BETTER!
			} catch (boost::program_options::unknown_option e) {
				cout << "Unrecognised option: " << e.what() << endl;
				cout << desc << endl;
				cout << "Exiting" << endl;
				exit(1); // NOT ELEGANT - DO BETTER!
			}

			po::notify(vm);

			if (vm.count("help")) {
			    cout << desc << endl;
			    exit(1); // NOT ELEGANT - DO BETTER!
			}

			if (vm.count("verbose")) {
				verbose = true;
			}

			// Process program options
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout << "PROGRAM OPTIONS" << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

			// Defaults are set above
			if (vm.count("source")) {
			    source = vm["source"].as<MCsource>();
			    if (source == signal) {
			    	doSignal = true;
			    } else {
			    	doSignal = false;
			    }
			} else {
			    cout << "MC source was not set. Defaulting to signal." << endl;
			}

			if (vm.count("swapMuRandomly")) {
			    swapMuRandomly = vm["swapMuRandomly"].as<bool>();
			} else {
			    cout << "Mu ordering not set. Defaulting to random." << endl;
			}

			if (vm.count("doHLT")) {
			    doHLT = vm["doHLT"].as<bool>();
			} else {
			    cout << "HLT requirement not set. Defaulting to using samples with HLT cuts." << endl;
			}

			if (vm.count("doRescale")) {
				doRescale = true;
				rescale = alpha;
			} else {
				cout << "No rescaling option set. Defulating to NOT rescaling to data." << endl;
			}

			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		} // end of constructor

		// Getters
		MCsource getSource() { return source; }
		bool getSignal() { return doSignal; }
		bool getQCDMu() { return doMu; }
		bool getMuOrdering() { return swapMuRandomly; }
		bool getHLT() { return doHLT; }
		int  getNEvents() { return nEvents; }
		bool getVerbose() { return verbose; }
		double getdR() { return dR; }
		double getRescale() { return rescale; }
		bool doRescaling() { return doRescale; }

		// This should really be in a separate .cc file...
		void printProgramOptions() {
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout << "PROGRAM OPTIONS" << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

			cout << "Doing " << source << " MC." << endl;

			if (swapMuRandomly) {
				cout << "Swapping mu 1<->2 randomly." << endl;
			} else {
				cout << "Mu 1,2 are pT ordered." << endl;
			}

			if (doHLT) {
				cout << "Using MC with HLT cuts already applied." << endl;
			} else {
				cout << "Using MC without any HLT cuts." << endl;
			}

			if (doRescale) {
				// cout << "Rescaling track eta & phi to match data. Factor: " << rescale << endl;
				cout << "Rescaling track eta & phi to match data. Doing quantile rescaling. " << endl;
			} else {
				cout << "Not rescaling track eta & phi to match data." << endl;
			}

			cout << "DeltaR(mu-mu) > " << dR << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}

};

/**
 * Add ROOT files to TChain from specified folder, from file_startNum to file_endNum, inclusive.
 * @param chain    TChain to add files to
 * @param folder   folder name
 * @param file     File name, of form myFile_whatever_
 * @param startNum Starting file number
 * @param endNum   Ending file number
 */
void addFilesFromFolder(TChain* chain, std::string folder, std::string file, int startNum, int endNum) {
	cout << "Adding " << folder+file+boost::lexical_cast<std::string>(startNum)+".root"
		 <<  " to " << folder+file+boost::lexical_cast<std::string>(endNum)+".root" << endl;
	for (int i = startNum; i <= endNum; i++) {
		chain->Add((folder+file+boost::lexical_cast<std::string>(i)+".root").c_str());
	}
}

/**
 * Add ROOT files to TChain from specified folder, from file_1 to file_endNum, inclusive.
 * @param chain    TChain to add files to
 * @param folder   folder name
 * @param file     File name, of form myFile_whatever_
 * @param endNum   Ending file number
 */
void addFilesFromFolder(TChain* chain, std::string folder, std::string file, int endNum) {
	addFilesFromFolder(chain, folder, file, 1, endNum);
}

/**
 * Add Delphes ntuples to TChain so we can process them in one go
 * @param chain    Pointer to TChain to add files to
 * @param pOpts    Pointer to ProgramOpts object that holds all user flags
 * (makes it easy to change avaiable flags without updating all programs individually!)
 */
void addInputFiles(TChain* chain, ProgramOpts* pOpts) {
	// Create chain of root trees
	// expect files to be named like myFile_i.root, where i = 1 -> nFiles

	MCsource source = pOpts->getSource();
	bool doMu       = pOpts->getQCDMu();
	bool doHLT      = pOpts->getHLT();

	////////////////
	// SIGNAL MC //
	////////////////
	if (source == signal) {
		if (doHLT) {
			cout << "Doing signal with HLT cuts" << endl;
			// folder = "Signal_1prong_500K_TauPythia_bare/";
			// file = "signal_1prong_500K_TauPythia_HLT_";
			// folder = "Signal_1prong_500K_bare/";
			// file = "signal_1prong_500K_HLT_";
			addFilesFromFolder(chain, "Signal_1prong_HLT_bare/", "Signal_HLT_", 60);
		} else {
			cout << "Doing signal without HLT cuts" << endl;
			addFilesFromFolder(chain, "Signal_1prong_500K_bare/", "signal_1prong_500K_NoHLT_", 20);
		}

	////////////////
	// QCDb MC	 //
	////////////////
	} else if (source == qcdb) {
		if (doMu) {
			if(doHLT) {
				cout << "Doing QCDb_mu with HLT cuts" << endl;
				addFilesFromFolder(chain, "QCDb_HLT_bare/", "QCDb_HLT_", 2000);
			} else{
				cout << "Doing QCDb_mu without HLT cuts" << endl;
				throw std::invalid_argument("DON'T USE QCDb_mu");
				addFilesFromFolder(chain, "QCDb_mu_pthatmin20_bare/", "QCDb_mu_pthatmin20_", 94);
			}
		} else {
			cout << "Doing QCDb" << endl;
			addFilesFromFolder(chain, "QCDb_cleanTk/", "QCDb_", 10);
		}

	//////////////
	// QCDc MC //
	//////////////
	} else if (source == qcdc) {
		if(doHLT) {
			cout << "Doing QCDc_mu with HLT cuts" << endl;
			addFilesFromFolder(chain, "QCDc_mu_pthatmin20_Mu17_Mu8_bare/", "QCDc_mu_pthatmin20_Mu17_Mu8_", 200);
		}

	//////////////////////
	// QCD q-g scatter //
	//////////////////////
	} else if (source == qcdscatter) {
		if(doHLT) {
			cout << "Doing QCDScatter with HLT cuts" << endl;
			addFilesFromFolder(chain, "QCDbcScatter_HLT_bare/", "QCDbcScatter_HLT_250_", 1, 200);
			addFilesFromFolder(chain, "QCDbcScatter_HLT_bare/", "QCDbcScatter_HLT_500_", 201, 800);
		}

	////////////////////
	// QCD all? OLD! //
	////////////////////
	} else if (source == qcdall) {
		if (doHLT) {
			cout << "Doing QCDAll with HLT cuts" << endl;
			// folder = "QCDAll_mu_pthatmin20_Mu17_Mu8_bare/";
			// file = "QCDAll_mu_pthatmin20_Mu17_Mu8_";
			addFilesFromFolder(chain, "QCDAll_NEW_mu_pthatmin20_Mu17_Mu8_bare/", "QCDAll_NEW_mu_pthatmin20_Mu17_Mu8_", 200);
		}

	////////////////
	// Test file //
	////////////////
	} else if (source == test) {
		cout << "Doing test files" << endl;
		addFilesFromFolder(chain, "QCDb_HLT_bare/", "QCDb_HLT_", 5);
	}

	// Auto-loop over ROOT files in folder using Boost::Filesystem

	// TODO
}

#endif
