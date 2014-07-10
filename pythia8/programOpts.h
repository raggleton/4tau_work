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
// -L$(HOME)/boost_1_55_0_install/lib -lboost_program_options 
// to complate make list (after -L$(HEPMCLOCATION)/lib -lHepMC \ )
#include <boost/lexical_cast.hpp>
// #include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using std::cout;
using std::endl;

// namespace fs = boost::filesystem;
namespace po = boost::program_options;

// Global enum for QCD process to be simulated
enum process { qcdb, qcdc, qcdscatter };
// Have to define << and >> ops to get enum to work with boost::lexical_cast 
// and program_options, and also so we can cout the enum easily
std::ostream& operator<<(std::ostream& out, const process value) {

    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value) {
        PROCESS_VAL(qcdb);     
        PROCESS_VAL(qcdc);
        PROCESS_VAL(qcdscatter);
    }
#undef PROCESS_VAL

    return out << s;
}

// Must be a way to improve this using static_cast
std::istream & operator>>(std::istream & in, process & value) {
  std::string token;
  if (in >> token) {
  	if (token == "qcdb")
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

class ProgramOpts
{
	private:
		bool outputEvent; 
		bool writeHLTToHEPMC; 
		bool writeNoHLTToHEPMC; 
		process userProcess;
		bool notMuOnly;
		int nEvents;
		std::string filename;
		bool verbose;
		bool disableScatterHook;
		int seed;

	public: 
		// constructor, parses input
		ProgramOpts(int argc, char* argv[]):
			outputEvent(false), 
			writeHLTToHEPMC(false), 
			writeNoHLTToHEPMC(false), 
			userProcess(qcdb),
			notMuOnly(false),
			nEvents(1),
			filename("testChangeMe"),
			verbose(false),
			disableScatterHook(false),
			seed(0)
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
				("process,p", po::value<process>(&userProcess),
					"Process to run: qcdb [default], qcdc, qcdscatter")
				("notMuOnly", 
					"Each event is kept regardless of # muons (default keeps hadronising event until 2+ muons)")
				("noScatterBC",
					"Make qcdscatter process use all flavours of q/qbar (defualt is b/c only)")
				("number,n", po::value<int>(&nEvents), 
					"Number of events to run over [default = 500]. If writeHLT enabled, counts # events passing HLT. Otherwise, counts # events with 2+ muons.")
				("name", po::value<std::string>(&filename),
					"Stem for output HepMC filenames (produces name_HLT.hepmc and name_NoHLT.hepmc)")
				("seed", po::value<int>(&seed),
					"Seed for random number generator. 0 = uses time. DON'T USE FOR BATCH SYSTEM. Get simultaneous start = same seed = pain. Use file numebr instead.")
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
			if (vm.count("notMuOnly")) {
			    notMuOnly = true;
			}
			if (vm.count("process")) {
			    userProcess = vm["process"].as<process>();
			}			
			if (vm.count("noScatterBC")) {
			    disableScatterHook = true;
			}

		} // end of constructor

		// Getters
		bool getOutputEvent() { return outputEvent; } 
		bool getWriteHLTToHEPMC() { return writeHLTToHEPMC; } 
		bool getWriteNoHLTToHEPMC() { return writeNoHLTToHEPMC; } 
		bool getNotMuOnly() { return notMuOnly; }
		process getProcess() { return userProcess; }
		int getNEvents() { return nEvents; }
		std::string getFilename() { return filename; }
		bool getVerbose() { return verbose; }
		bool getScatterHook() { return disableScatterHook; }
		int getSeed() { return seed; }

		// This should really be in a separate .cc file...
		void printProgramOptions() {
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout << "PROGRAM OPTIONS" << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			
			cout << "Process: " << userProcess << endl;
			if (outputEvent)
				cout << "Outputting first event" << endl;
			if (writeHLTToHEPMC)
				cout << "Writing events that pass HLT mu cuts to hepmc file" << endl;
			if (writeNoHLTToHEPMC)
				cout << "Writing events that have >= 2 muons to hepmc file" << endl;
			if (notMuOnly)
				cout << "Don't care if 2+ muons" << endl;
			if (disableScatterHook)
				cout << "q-g scatter: scatter process use all flavours of q/qbar" << endl;
			cout << "Doing " << nEvents << " events" << endl;
			cout << "Writing to " << filename <<"(_HLT|_NoHLT).hepmc" << endl;
			cout << "Random seed: " << seed << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}

};

#endif
