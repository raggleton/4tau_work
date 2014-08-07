// Based off main15.cc in examples folder
// 
// This creates QCD b samples, which can optionally decay to muons

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"   
#include "HepMC/IO_GenEvent.h"
#include "myHooks.h"
#include "programOpts.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace Pythia8;
using namespace boost::algorithm;

std::string getCurrentTime() {
    // current date/time based on current system
  time_t now = time(0);
  // convert now to string form
  char* dt = ctime(&now);
  std::string str1 = std::string(dt);
  trim(str1);
  return str1;
}
 
int main(int argc, char* argv[]) {

  ProgramOpts pOpts(argc, argv);
  pOpts.printProgramOptions();

  bool outputEvent       = pOpts.getOutputEvent(); // output entire event listing to STDOUT (long!), for debugging only
  bool writeHLTToHEPMC   = pOpts.getWriteHLTToHEPMC(); // output to HEPMC events passing HLT
  bool writeNoHLTToHEPMC = pOpts.getWriteNoHLTToHEPMC(); // output to HEPMC events without any HLT cuts
  bool notMuOnly         = pOpts.getNotMuOnly(); // Ture is the user doesn't want events generated that are guaranteed to have 2+ muons
  bool DEBUG             = pOpts.getVerbose();
  int seed               = pOpts.getSeed();
  process userProcess    = pOpts.getProcess(); // qcdb, qcdc, qcdscatter
  std::string filename   = pOpts.getFilename(); // HEPMC filename stem to be used

  if (writeHLTToHEPMC) cout << "Writing HLT events to HEPMC" << endl;
  if (writeNoHLTToHEPMC) cout << "Writing all events to HEPMC" << endl;

  // Interface for conversion from Pythia8::Event to HepMC event. 
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  // Do one for with HLT cuts, one wihtout HLT cuts
  std::string noHLTfile = filename+"_NoHLT.hepmc";
  std::string HLTfile = filename+"_HLT.hepmc";
  HepMC::IO_GenEvent ascii_io_NoHLT(noHLTfile, std::ios::out);
  HepMC::IO_GenEvent ascii_io_HLT(HLTfile, std::ios::out);

  // Text file to write progress - handy for monitoring during PBS jobs
  ofstream myfile;
  myfile.open((filename+"_progress.txt").c_str());

  cout << "Outputting to " << noHLTfile << " and " << HLTfile << endl;
  cout << "Outputting iEvent status to " << filename << ".txt" << endl;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;  
  Event savedEvent; // used for holding event when doing repeat hadronisation
  
  // Do random seed
  // A value 0 gives a random seed based on the time
  pythia.readString("Random:setSeed = on");
  std::string seedStr = "Random:seed = "+boost::lexical_cast<std::string>(seed);
  pythia.readString(seedStr);
  
  ////////////////////////
  // Number of events.  //
  ////////////////////////
  // For HLT. NoHLT has about 60X HLT amount (320K NoHLT evt for 5K HLT evnt)
  // Warning, 5K events ~900MB hepmc file and takes ~5 min.
  // Warning, 50K events ~9GB hepmc file and takes ~40 min.
  // pythia.readString("Main:numberOfEvents = 100");
  int nEvent = pOpts.getNEvents();
  pythia.readString("Next:numberShowEvent = 00");

  ///////////////////////////
  // Setup for pp at 8 TeV //
  ///////////////////////////
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 8000");

  // Making it like PythiaUEZ2starSettings
  pythia.readString("PDF:pSet = 8"); // CTEQ6L1

  ////////////////////////////////////////////////////
  // Simulate production above given pTmin scale. //
  ////////////////////////////////////////////////////
  // Warning: these processes do not catch all possible production modes.
  // You would need to use HardQCD:all or even SoftQCD:nonDiffractive for that.
  
  // g-q scatter:
  if(userProcess == qcdscatter) {
    pythia.readString("HardQCD:nQuarkNew = 5");    
    pythia.readString("HardQCD:qg2qqqbarDiff = on");    
    pythia.readString("HardQCD:qg2qqqbarSame= on");    
    QGScatterHook* scatterHook = new QGScatterHook(DEBUG); // should smart pointer this
    if (!pOpts.getScatterHook()){
      pythia.setUserHooksPtr(scatterHook);
    }
  }
  
  // gg/qqbar - > ccbar:
  if (userProcess == qcdc) {
    pythia.readString("HardQCD:hardccbar = on");    
  }
  
  // gg/qqbar - > bbbar:
  if (userProcess == qcdb) {
    pythia.readString("HardQCD:hardbbbar = on");    
  }

  // Extra ttbar contribution
  // pythia.readString("Top:gg2ttbar = on");
  // pythia.readString("Top:qqbar2ttbar = on");
  // Make sure t->Wb only
  // pythia.readString("6: onMode = off");
  // pythia.readString("6: onIfAny = 5");

  pythia.readString("PhaseSpace:pTHatMin = 20.");  

  // pythia.readString("ProcessLevel:all = off");   
  // pythia.readString("PartonLevel:all = off");   
  
  // Do repeated hadronisation to get 2+ muons out of each event
  if (!pOpts.getNotMuOnly()) {
    pythia.readString("HadronLevel:all = off");
  }

  // Initialize 
  pythia.init();

  // Some basic histograms
  Hist nMuInEvent("number of muons in an event", 10, -0.5, 9.5); 
  Hist nMuInEventHLT("number of muons in an event (HLT)", 10, -0.5, 9.5); 
  Hist nMuInEvent2plus("number of muons in an event (2+ selection)", 10, -0.5, 9.5);

  Hist muPt("pT muons in an event (HLT)", 40, 0.0, 40.0); 
  Hist muPtNoHLT("pT muons in an event (2+ muons)", 40, 0.0, 40.0); 
  
  Hist nRepeats("number of repeats required to get 2+ muons", 50, 0.5, 100.5); 
  int nWithSSMuPair = 0;

  cout << "Starting at " << getCurrentTime() << endl;
  myfile << "Starting at " << getCurrentTime() << endl;

  ///////////////////////
  // Begin event loop. //
  ///////////////////////
  int iEvent      = 0; // counts events passing HLT, NOT NoHLT
  int lastiEvent  = 0; // for outputting to screen every outputEvery events
  int outputEvery = 5; 
  bool printEvent = true;
  int counter = 0;
  while(iEvent < nEvent) {
    bool wantedHLT = false;
    bool wantedNoHLT = false;
    counter++;
    if ((iEvent % outputEvery == 0) && (iEvent!= lastiEvent)){
      lastiEvent = iEvent;
      cout << "iEvent: " << iEvent << " - " << getCurrentTime() << endl;
      myfile << "iEvent: " << iEvent << " - " << getCurrentTime() << endl;
    }
   
    // Generate event. Skip it if error.
    if (!pythia.next()) continue;

    ///////////////////////////////////////////////////////////////////
    // Now do all your selections to determine if we keep the event: //
    ///////////////////////////////////////////////////////////////////

    savedEvent = event;

    // Look for muons among decay products (also from charm/tau/...).
    int nMuNeg(0), nMuPos(0);
    std::vector<double> muPtVec;
    
    
    for (int i = 0; i < event.size(); ++i) {
      int id = event[i].id();  
      int status = event[i].status();
      if (id ==  13 && status > 0){ 
        nMuNeg++;
        muPtVec.push_back(event[i].pT());
      }
      if (id == -13 && status > 0) {
        nMuPos++;
        muPtVec.push_back(event[i].pT());
      }
    }
    
    int iRepeat = 0; // count repeated decays of same event
    // Unless user wants any # muons in their events, 
    // we continually hadronise until we get 2+ muons
    if (!pOpts.getNotMuOnly()){
      int nAttempts = 0;
      while(muPtVec.size() < 2 && nAttempts <100) {
        nAttempts++;
        muPtVec.clear();
        nMuNeg = 0;
        nMuPos = 0;

        if (iRepeat > 0) event = savedEvent;

        iRepeat++;
        if (!pythia.forceHadronLevel(false)) continue;


        for (int i = 0; i < event.size(); ++i) {
          int id = event[i].id();  
          int status = event[i].status();
          if (id ==  13 && status > 0){ 
            nMuNeg++;
            muPtVec.push_back(event[i].pT());
          }
          if (id == -13 && status > 0) {
            nMuPos++;
            muPtVec.push_back(event[i].pT());
          }
        }
        nMuInEvent.fill(muPtVec.size());
      }

      // if (DEBUG) event.list();
      nRepeats.fill(iRepeat);
      counter += iRepeat;
    }

    // Check whether SS pair(s) present.
    if ((nMuNeg  > 1) || (nMuPos > 1)) {
      ++nWithSSMuPair;
    }

    nMuInEvent2plus.fill(muPtVec.size());

    if (!writeHLTToHEPMC) iEvent++; // Count evt if not doing HLT mode. 
    
    wantedNoHLT = true;
    
    for (unsigned a = 0; a < muPtVec.size(); a++){
        muPtNoHLT.fill(muPtVec.at(a));
    }
    
    // Do HLT cuts
    if (writeHLTToHEPMC){
      
      if (muPtVec.size() != 0){
        // order mu pt vector
        std::sort(muPtVec.begin(),muPtVec.end(), std::greater<int>());

        // Emulate HLT - HLT_Mu17_Mu8
        if (muPtVec[0] > 17 && muPtVec[1] > 8 ){
          // if it gets to here, then we're happy with the event
          iEvent++;
          for (unsigned a = 0; a < muPtVec.size(); a++){
            muPt.fill(muPtVec.at(a));
          }
          nMuInEventHLT.fill(nMuPos + nMuNeg);
          wantedHLT = true;
        }
      }
    }

    // Output the first HLT event to screen
    if (outputEvent && iEvent == 1 && printEvent){
      event.list();
      printEvent = false;
    }

    // Write out events that pass HLT cuts
    if (wantedHLT && writeHLTToHEPMC){
      // Construct new empty HepMC event and fill it.
      // Units will be as chosen for HepMC build, but can be changed
      // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)  
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      ToHepMC.fill_next_event( pythia, hepmcevt );

      // Write the HepMC event to file. Done with it.
      ascii_io_HLT << hepmcevt;
      delete hepmcevt;
    }

    // Write out events that have 2+ muons, regardless of whether they pass HLT cuts
    if (wantedNoHLT && writeNoHLTToHEPMC){
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      ToHepMC.fill_next_event( pythia, hepmcevt );

      // Write the HepMC event to file. Done with it.
      ascii_io_NoHLT << hepmcevt;
      delete hepmcevt;
    }
  } //end of events loop

  // Statistics. Histograms. 
  pythia.stat();
  cout << muPt << nMuInEvent << nMuInEvent2plus << nMuInEventHLT << muPtNoHLT << nRepeats << endl;
  cout << "Number of events with pair of SS muon: " << nWithSSMuPair << endl;
  cout << counter << endl;
  // Done.
  // if(userProcess == qcdscatter) {
    // delete scatterHook;
  // }
  myfile.close();
  return 0;
}
