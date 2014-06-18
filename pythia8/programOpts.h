#ifndef PROGRAMOPTS_H
#define PROGRAMOPTS_H

// #include <vector>
// #include <algorithm>
#include <string>
#include <sstream>
#include <iostream>

// BOOST headers
// Need to add 
// -I $(HOME)/boost_1_55_0 -I $(HOME)/boost_1_55_0_install/include 
// to CXXFLAGS in Delphes/Makefile
// and
// -L/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib -lboost_program_options -lboost_filesystem
// to complate make list (after -L$(HEPMCLOCATION)/lib -lHepMC \ )
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using std::cout;
using std::endl;

namespace fs = boost::filesystem;
namespace po = boost::program_options;

class ProgramOpts
{
	private:
		bool outputEvent; 
		bool writeHLTToHEPMC; 
		bool writeNoHLTToHEPMC; 
		bool muOnly;
		bool tauToMuOnly;
		int nEvents;
		std::string filename;
		bool verbose;

	public: 
		// constructor, parses input
		ProgramOpts(int argc, char* argv[]):
			outputEvent(false), 
			writeHLTToHEPMC(false), 
			writeNoHLTToHEPMC(false), 
			muOnly(false),
			tauToMuOnly(false),
			nEvents(500),
			filename("testChangeMe"),
			verbose(false)
		{
			po::options_description desc("\nAllowed options");
			desc.add_options()
				("help,h", "Produce help message")
				("outputEvent", 
					"Outputs complete event listing of first event passing HLT (if writeHLT enabled) or with >= 2 muons")
				("writeHLT",  
					"write events passing HLT mu cuts to file")
				("writeNoHLT",
					"write events with >= 2 muons to file")
				("muOnly", 
					"Allow b/c hadrons to only decay to final state which contains a muon")
				("tauToMuOnly", 
					"Enforce Taus from b/c hadrons decays to decay to muons")
				("number,n", po::value<int>(&nEvents), 
					"Number of events to run over [default = 500]")
				("name", po::value<std::string>(&filename),
					"Stem for output HepMC filenames (produces name_HLT.hepmc and name_NoHLT.hepmc)")
				("verbose,v", "Output debugging statements")
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

			
			// Defaults are set above
			if (vm.count("outputEvent")) {
			    outputEvent = true;
			}
			if (vm.count("writeHLT")) {
			    writeHLTToHEPMC = true;
			}
			if (vm.count("writeNoHLT")) {
			    writeNoHLTToHEPMC = true;
			}
			if (vm.count("muOnly")) {
			    muOnly = true;
			}
			if (vm.count("tauToMuOnly")) {
			    tauToMuOnly = true;
			}

		} // end of constructor

		// Getters
		bool getOutputEvent() { return outputEvent; } 
		bool getWriteHLTToHEPMC() { return writeHLTToHEPMC; } 
		bool getWriteNoHLTToHEPMC() { return writeNoHLTToHEPMC; } 
		bool getMuOnly() { return muOnly; }
		bool getTauToMuOnly() { return tauToMuOnly; }
		int getNEvents() { return nEvents; }
		std::string getFilename() { return filename; }
		bool getVerbose() { return verbose; }

		// This should really be in a separate .cc file...
		void printProgramOptions() {
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout << "PROGRAM OPTIONS" << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			
			if (outputEvent)
				cout << "Outputting first event" << endl;
			if (writeHLTToHEPMC)
				cout << "Writing events that pass HLT mu cuts to hepmc file" << endl;
			if (writeNoHLTToHEPMC)
				cout << "Writing events that have >= 2 muons to hepmc file" << endl;
			if (muOnly)
				cout << "Forcing b/c hadrons to decay to semi-muonic final state" << endl;
			if (tauToMuOnly)
				cout << "Force taus from b/c hadrons to decay to muons" << endl;
			cout << "Doing " << nEvents << " events" << endl;
			cout << "Writing to " << filename <<"(_HLT|_NoHLT).hepmc" << endl;
			
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}

};

#endif
