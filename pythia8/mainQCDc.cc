// Based off main15.cc in examples folder
// 
// This creates QCD c samples, which can optionally decay to muons

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"   
#include "HepMC/IO_GenEvent.h"

using namespace Pythia8;
 
int main(int argc, char* argv[]) {

  bool outputEvent       = false; // output entire event listing to STDOUT (long!), for debugging only
  bool writeHLTToHEPMC   = true; // output to HEPMC events passing HLT
  bool writeNoHLTToHEPMC = false; // output to HEPMC events without any HLT cuts
  bool muOnly            = true; // Only allow b hadrons to decay to muons or taus
  bool tauToMuOnly       = true; // Only allow those taus from b hadrons to decay to muons 

  // Check that correct number of command-line arguments
  // Unfortunately required even if writeToHEPMC = false
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments. \n "
         <<  "You are expected to provide one output file name eg myQCDb \n"
         << " Program stopped! " << endl;
    return 1;
  }

  if (writeHLTToHEPMC) cout << "Writing HLT events to HEPMC" << endl;
  if (writeNoHLTToHEPMC) cout << "Writing all events to HEPMC" << endl;

  // Interface for conversion from Pythia8::Event to HepMC event. 
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  // Do one for with HLT cuts, one wihtout HLT cuts
  std::string noHLTfile = std::string(argv[1])+"_NoHLT.hepmc";
  std::string HLTfile = std::string(argv[1])+"_HLT.hepmc";
  HepMC::IO_GenEvent ascii_io_NoHLT(noHLTfile, std::ios::out);
  HepMC::IO_GenEvent ascii_io_HLT(HLTfile, std::ios::out);

  cout << "Outputting to " << noHLTfile << " and " << HLTfile << endl;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;
  
  // Do random seed
  //  a value 0 gives a random seed based on the time
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  
  // Number of events. For HLT. NoHLT has about 60X HLT amount (320K NoHLT evt for 5K HLT evnt)
  // Warning, 5K events ~900MB hepmc file and takes ~5 min.
  // Warning, 50K events ~9GB hepmc file and takes ~40 min.
  pythia.readString("Main:numberOfEvents = 5000");
  int nEvent = pythia.mode("Main:numberOfEvents");
  pythia.readString("Next:numberShowEvent = 00");
  // pythia.readString("Next:numberShowProcess = 100");
  

  // Setup for pp at 8 TeV
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 8000");

  // Making it like PythiaUEZ2starSettings
  pythia.readString("PDF:pSet = 8"); // CTEQ6L1
  // pythia.readString("");

  // Simulate ccbar production above given pTmin scale.
  // Warning: these processes do not catch all possible production modes.
  // You would need to use HardQCD:all or even SoftQCD:nonDiffractive for that.
  pythia.readString("HardQCD:gg2ccbar = on");    
  pythia.readString("HardQCD:qqbar2ccbar = on");    
  // Make sure t->Wb only
  pythia.readString("PhaseSpace:pTHatMin = 20.");  
  // pythia.readString("HadronLevel:all = on");
    // pythia.readString("ProcessLevel:all = off");   
  // pythia.readString("PartonLevel:all = off");   
  // pythia.readString("HadronLevel:all = off");   

  // All C hadrons from PDG
  std::vector<int> cCodes{ 411,421,10411,10421,413,423,10413,10423,20413,20423,415,425,431,10431,433,
    10433,20433, 435,441,10441,100441,443,10443,20443,100443,30443,9000443,9010443,9020443,445,100445,
    4122,4222,4212,4112,4224,4214,4114,4232,4132,4322,4312,4324,4314,4332,4334,4412,4422,4414,4424,4432,4434,4444 };
  int nCodes = cCodes.size();

  if (muOnly){
    // For D hadrons to decay weakly to muons or taus
    for (int iC = 0; iC < nCodes; ++iC) {
      // Check PDGID is in PYTHIA
      if(! pythia.particleData.isParticle(cCodes[iC])) continue;
      
      // Get particle name.
      // If excited state, then just skip it, we want it to decay to less excited states.
      if (pythia.particleData.name(cCodes[iC]).find("*") != std::string::npos) continue;
      
      // Turn off all decay modes first
      std::stringstream sstm;
      sstm << cCodes[iC] << ":onMode = off";
      std::string command = sstm.str();
      pythia.readString(command);
      sstm.str("");
      
      // Now just turn on tau or mu ones
      sstm << cCodes[iC] << ":onIfAny = 13 15";
      command = sstm.str();
      pythia.readString(command);
    }
    // pythia.particleData.list(511);
  }

  // For tau, turn off decays. We want any tau from B hadrons to decay to muons,
  // but all other taus can decay however they want.
  if(tauToMuOnly)
    pythia.readString("15:onMode = off");
  
  std::vector<int> tausFromC;
  std::vector<int> tausNotFromC;

  // Initialize 
  pythia.init();

  // Some basic historgrams
  Hist nMuInEvent("number of muons in an event (HLT)", 10, -0.5, 9.5); 
  Hist muPt("pT muons in an event (HLT)", 40, 0.0, 40.0); 
  Hist muPtNoHLT("pT muons in an event (NoHLT)", 40, 0.0, 40.0); 
  int nWithPair = 0;


  // Begin event loop.
  int iEvent = 0; // tracks events passing HLT, NOT NoHLT
  int lastiEvent = 0;
  while(iEvent < nEvent) {
    bool wantedHLT = false;
    bool wantedNoHLT = false;

    if ((iEvent % 50 == 0) && (iEvent!= lastiEvent)){
      lastiEvent = iEvent;
      cout << "iEvent: " << iEvent << endl;
    }

    if (tauToMuOnly){
      // have to turn off tau decays for each event
      pythia.readString("15:onMode = off"); 
      // reset tau vectors to empty
      tausFromC.resize(0);
      tausNotFromC.resize(0);
    }

    // Generate event. Skip it if error.
    if (!pythia.next()) continue;
    
    if (tauToMuOnly) {
      // begin particles-in-event loop
      // Let's find all the taus, and store their positions.
      // We want to store those from B hadrons differently to those not from B hadrons
      // and then decay them later
      for (int i = 0; i < event.size(); ++i) {
        
        // int stat = event[i].statusAbs();
        if (event[i].idAbs() == 15 ){ //&& stat == 91){
          // Mark as decayed 
          event[i].statusNeg();

          // Loop through and see if mother of tau is in C hadron list
          // Must be a better way?
          bool daughterOfC = false;
          // motherN() gets the particle # not the PDGID number
          for (int iC = 0; iC < nCodes; ++iC) {
            if (event[event[i].mother1()].idAbs() == cCodes[iC] || event[event[i].mother2()].idAbs() == cCodes[iC]){
              daughterOfC = true;
              event[i].statusNeg();
              break;
            }
          }

          if (daughterOfC)
            tausFromC.push_back(i);    
          else
            tausNotFromC.push_back(i);    

        } 
      } //end loop over particles

      // Now we do decays in 2 goes.
      // First do all the taus from Bs
      pythia.readString("15:onIfAny = 13");
      for(unsigned a =0; a < tausFromC.size(); a++){
        event[tausFromC[a]].statusPos();
      }
      if (!pythia.moreDecays()) continue; // Go ahead and decay the undecayed taus
      
      // Now do all other taus 
      // Turn on all modes of decay
      pythia.readString("15:onMode = on");
      for(unsigned a =0; a < tausNotFromC.size(); a++){
        event[tausNotFromC[a]].statusPos(); // Mark the taus not from B as needing decaying
      }
      if (!pythia.moreDecays()) continue;
    } // end of if(tauToMuOnly)
    
    ///////////////////////////////////////////////////////////////////////////////
    // Now do all your selection requirements to determine if we keep the event: //
    ///////////////////////////////////////////////////////////////////////////////
    
    // Look for muons among decay products (also from charm/tau/...).
    int nMuNeg(0), nMuPos(0);
    std::vector<double> muPtVec;

    for (int i = 0; i < event.size(); ++i) {
      int id = event[i].id();  
      if (id ==  13){ 
        nMuNeg++;
        muPtVec.push_back(event[i].pT());
      }
      if (id == -13) {
        nMuPos++;
        muPtVec.push_back(event[i].pT());
      }
    }
    
    // Check whether SS pair(s) present.
    if ((nMuNeg  > 1) || (nMuPos > 1)) {
      ++nWithPair;
    }

    if (nMuPos+nMuNeg < 2) continue; // Skip if there's only 1 muon
    
    if (!writeHLTToHEPMC) iEvent++; // Count evt if not doing HLT mode. Need to be careful when testing this without either...
    wantedNoHLT = true;
    
    for (unsigned a = 0; a < muPtVec.size(); a++){
        muPtNoHLT.fill(muPtVec.at(a));
    }
    
    // Do HLT cuts
    if (writeHLTToHEPMC){
      // order mu pt vector
      std::sort(muPtVec.begin(),muPtVec.end(), std::greater<int>());

      // Emulate HLT - HLT_Mu17_Mu8
      if (muPtVec[0] > 17 && muPtVec[1] > 8 ){
        // if it gets to here, then we're happy with the event
        iEvent++;
        for (unsigned a = 0; a < muPtVec.size(); a++){
          muPt.fill(muPtVec.at(a));
        }
        nMuInEvent.fill(nMuPos + nMuNeg);
        wantedHLT = true;
      }
    }

    // Output the event to screen
    if (outputEvent)
      event.list();

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
      // Construct new empty HepMC event and fill it.
      // Units will be as chosen for HepMC build, but can be changed
      // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)  
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      ToHepMC.fill_next_event( pythia, hepmcevt );

      // Write the HepMC event to file. Done with it.
      ascii_io_NoHLT << hepmcevt;
      delete hepmcevt;
    }
  } //end of events loop

  // Statistics. Histograms. 
  pythia.stat();
  cout << muPt << nMuInEvent << muPtNoHLT << endl;
  cout << "Number of events with pair, & passing HLT: " << nWithPair << endl;

  // Done. 
  return 0;
}
