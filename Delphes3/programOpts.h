#ifndef PROGRAMOPTS_H
#define PROGRAMOPTS_H

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>

// ROOT headers
// #include "TStyle.h"

// BOOST headers
// Need to add 
// -I $(HOME)/boost_1_55_0 -I $(HOME)/boost_1_55_0_install/include 
// to CXXFLAGS in Delphes/Makefile
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

// Global enum for the MC source
// Note convention of all lower case
enum MCsource { signal, qcdb, qcdc, qcdscatter };

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
		MCsource source; // do signal or qcd(b)(c)(scatter)
		bool doSignal; // do signal, or not signal (some parts are signal only)
		bool doMu; // for QCDb - either inclusive decays or mu only decays - DEPRECIATED
		bool swapMuRandomly; // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
		bool doHLT; // whether to use MC that has HLT cuts already applied.
		int  nEvents; // how many events to run over, -1 = all events

	public: 
		// constructor, parses input
		ProgramOpts(int argc, char* argv[]):
			// some sensible defaults
			source(signal),
			doSignal(true),
			doMu(true),
			swapMuRandomly(true),
			doHLT(true),
			nEvents(-1)
		{
			po::options_description desc("Allowed options");
			desc.add_options()
				("help", "produce help message")
				("source", po::value<MCsource>(&source), 
					"Process to run: signal [default], qcdb, qcdc, qcdscatter")
				("swapMuRandomly", po::value<bool>(&swapMuRandomly), 
					"TRUE [default] - mu 1,2 randomly assigned, FALSE - mu 1,2 pT ordered")
				("doHLT", po::value<bool>(&doHLT), 
					"TRUE [default] - use samples with HLT_Mu17_Mu8 during generation, FALSE - no HLT cuts")
				("nnumber,N", po::value<int>(&nEvents), 
					"Number of events to run over. -1 for all [defualt]")
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
			
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		} // end of constructor

		// Getters
		MCsource getSource() { return source; }
		bool getSignal() { return doSignal; }
		bool getQCDMu() { return doMu; }
		bool getMuOrdering() { return swapMuRandomly; }
		bool getHLT() { return doHLT; }
		int  getNEvents() { return nEvents; }
		
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
			
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}

};

/**
 * Add Delphes ntuples to TChain so we can process them in one go
 * @param chain    Pointer to TChain to add files to
 * @param pOpts    Pointer to ProgramOpts object that holds all user flags 
 * (makes it easy to change avaiable flags without updating all programs individually!)
 */
void addInputFiles(TChain* chain, ProgramOpts* pOpts) {
	// Create chain of root trees
	int nFiles = 0; // number of files to be added
	std::string folder = ""; // folder & file stem 
	std::string file = ""; // folder & file stem 
	// expect files to be named like myFile_i.root, where i = 1 -> nFiles
	
	MCsource source = pOpts->getSource();
	bool doMu       = pOpts->getQCDMu();
	bool doHLT      = pOpts->getHLT();

	if (source == signal) {
		if (doHLT) {
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_1_HLT_bare.root");
			cout << "Doing signal with HLT cuts" << endl;
			// folder = "Signal_1prong_500K_TauPythia_bare/";
			// file = "signal_1prong_500K_TauPythia_HLT_";
			folder = "Signal_1prong_500K_bare/";
			file = "signal_1prong_500K_HLT_";
			nFiles = 20;
			// nFiles = 1;
		} else { 
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_1_NoHLT_bare.root");
			cout << "Doing signal without HLT cuts" << endl;
			folder = "Signal_1prong_500K_bare/";
			file = "signal_1prong_500K_NoHLT_";
			nFiles = 20;
		}
	} else if (source == qcdb) {
		if (doMu) {
			if(doHLT) {
				cout << "Doing QCDb_mu with HLT cuts" << endl;
				folder = "QCDb_mu_pthatmin20_Mu17_Mu8_bare/";
				file = "QCDb_mu_pthatmin20_Mu17_Mu8_";
				nFiles = 350;
			} else{
				cout << "Doing QCDb_mu without HLT cuts" << endl;
				folder = "QCDb_mu_pthatmin20_bare/";
				file = "QCDb_mu_pthatmin20_";
				// chain->Add("QCDb_mu_pthatmin20_bare/QCDb_mu_pthatmin20_92.root");
				nFiles = 94;
			}
		} else {
			cout << "Doing QCDb" << endl;
			folder = "QCDb_cleanTk/";
			file = "QCDb_";
			nFiles = 10;
		}
	} else if (source == qcdc) {
		if(doHLT) {
			cout << "Doing QCDc_mu with HLT cuts" << endl;
			folder = "QCDc_mu_pthatmin20_Mu17_Mu8_bare/";
			file = "QCDc_mu_pthatmin20_Mu17_Mu8_";
			nFiles = 200;
		}
	} else if (source == qcdscatter) {
		if(doHLT) {
			cout << "Doing QCDScatter with HLT cuts" << endl;
			folder = "QCDScatter_mu_pthatmin20_Mu17_Mu8_bare/";
			file = "QCDScatter_mu_pthatmin20_Mu17_Mu8_";
			nFiles = 11;
		}
	}

	// Auto-loop over ROOT files in folder using Boost::Filesystem
	
	// TODO


	// For manually looping over files in a folder from 1 to nFiles (inclusive)	
	for (int i = 1; i <= nFiles; i ++) {
		cout << "Adding " << 
			folder+file+boost::lexical_cast<std::string>(i)+".root" << endl;
		chain->Add((folder+file+boost::lexical_cast<std::string>(i)+".root").c_str());
	}
}

#endif