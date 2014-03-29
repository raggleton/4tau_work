#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TArc.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

//------------------------------------------------------------------------------

// Instructions: in your Delphes folder do
// ./configure
// make (may need a make clean first)
// ./testScript_cleanTk (uses the filename of this cpp as the program name)

//------------------------------------------------------------------------------

// Here you can put your analysis macro   

#include "testScript_cleanTk.C"

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
	char *appName = "testScript_cleanTk";

	// if(argc != 2)
	// {
	//   cout << " Usage: " << appName << " input_file" << endl;
	//   cout << " input_file - input file in ROOT format ('Delphes' tree)," << endl;
	//   return 1;
	// }

	gROOT->SetBatch();

	int appargc = 1;
	char *appargv[] = {appName};
	TApplication app(appName, &appargc, appargv);

	TString inputFile(argv[1]);

//------------------------------------------------------------------------------

// Here you call your macro's main function 

	testScript_cleanTk();

//------------------------------------------------------------------------------

}


