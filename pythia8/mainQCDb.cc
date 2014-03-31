// Based off main15.cc in examples folder
// 
// This creates QCD b samples, which can optionally decay to muons

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"   
#include "HepMC/IO_GenEvent.h"

using namespace Pythia8;
 
int main(int argc, char* argv[]) {

  bool outputEvent  = false; // output entire event listing to STDOUT (long!), for debugging only
  bool writeToHEPMC = true; // output to HEPMC
  bool muOnly       = true; // Only allow b hadrons to decay to muons or taus
  bool tauToMuOnly  = true; // Only allow those taus from b hadrons to decay to muons 

  // Check that correct number of command-line arguments
  // Unfortunately required even if writeToHEPMC = false
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments. \n "
         <<  "You are expected to provide one output file name. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  cout << "Outputting to " << argv[1] << endl;

  // Interface for conversion from Pythia8::Event to HepMC event. 
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[1], std::ios::out);

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;
  
  // Do random seed
  //  a value 0 gives a random seed based on the time
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  
  // Number of events. 
  // Warning, 5K events ~900MB hepmc file and takes ~5 min.
  // Warning, 50K events ~9GB hepmc file and takes ~40 min.
  pythia.readString("Main:numberOfEvents = 50000");
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

  // Simulate b production above given pTmin scale.
  // Warning: these processes do not catch all possible production modes.
  // You would need to use HardQCD:all or even SoftQCD:nonDiffractive for that.
  pythia.readString("HardQCD:gg2bbbar = on");    
  pythia.readString("HardQCD:qqbar2bbbar = on");    
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");
  // Make sure t->Wb only
  pythia.readString("6: onMode = off");
  pythia.readString("6: onIfAny = 5");
  pythia.readString("PhaseSpace:pTHatMin = 20.");  
  // pythia.readString("HadronLevel:all = on");
    // pythia.readString("ProcessLevel:all = off");   
  // pythia.readString("PartonLevel:all = off");   
  // pythia.readString("HadronLevel:all = off");   

  // All B hadrons from PDG
  int bCodes[90] = {511,521,10511,10521,513,523,10513,10523,20513,20523,515,525,531,10531,533,10533,
    20533,535,541,10541,543,10543,20543,545,51,10551,100551,110551,200551,210551,553,10553,20553,30553,
    100553,110553,120553,130553,200553,210553,220553,300553,10860,9000553,11020,9010553,555,10555,20555,
    100555,110555,120555,200555,557,100557,5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,5314,
    5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,5514,5524,5532,5534,
    5542,5544,5554};
  int nCodes = 90;

  if (muOnly){
    // For B hadrons to decay weakly to muons or taus
    for (int iC = 0; iC < nCodes; ++iC) {
      // Check PDGID is in PYTHIA
      if(! pythia.particleData.isParticle(bCodes[iC])) continue;
      
      // Get particle name.
      // If excited state, then just skip it, we want it to decay to less excited states.
      if (pythia.particleData.name(bCodes[iC]).find("*") != std::string::npos) continue;
      
      // Turn off all decay modes first
      std::stringstream sstm;
      sstm << bCodes[iC] << ":onMode = off";
      std::string command = sstm.str();
      pythia.readString(command);
      sstm.str("");
      
      // Now just turn on tau or mu ones
      sstm << bCodes[iC] << ":onIfAny = 13 15";
      command = sstm.str();
      pythia.readString(command);
    }
    // pythia.particleData.list(511);
  }

  // For tau, turn off decays. We want any tau from B hadrons to decay to muons,
  // but all other taus can decay however they want.
  if(tauToMuOnly)
    pythia.readString("15:onMode = off");
  
  std::vector<int> tausFromB;
  std::vector<int> tausNotFromB;

  // Initialize 
  pythia.init();

  // Some basic historgrams
  Hist nMuInEvent("number of muons in an event", 10, -0.5, 9.5); 
  Hist muPt("pT muons in an event", 40, 0.0, 20.0); 
  int nWithPair = 0;


  // Begin event loop.
  // for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
  int iEvent = 0;
  int lastiEvent = 0;
  while(iEvent < nEvent) {
    bool wanted = true;

    if ((iEvent % 1000 == 0) && (iEvent!= lastiEvent)){
      lastiEvent = iEvent;
      cout << "iEvent: " << iEvent << endl;
    }

    if (tauToMuOnly){
      // have to turn off tau decays for each event
      pythia.readString("15:onMode = off"); 
      // reset tau vectors to empty
      tausFromB.resize(0);
      tausNotFromB.resize(0);
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

          // Loop through and see if mother of tau is in B hadron list
          // Must be a better way?
          bool daughterOfB = false;
          // motherN() gets the particle # not the PDGID number
          for (int iC = 0; iC < nCodes; ++iC) {
            if (event[event[i].mother1()].idAbs() == bCodes[iC] || event[event[i].mother2()].idAbs() == bCodes[iC]){
              daughterOfB = true;
              event[i].statusNeg();
              break;
            }
          }

          if (daughterOfB)
            tausFromB.push_back(i);    
          else
            tausNotFromB.push_back(i);    

        } 
      } //end loop over particles

      // Now we do decays in 2 goes.
      // First do all the taus from Bs
      pythia.readString("15:onIfAny = 13");
      for(unsigned a =0; a < tausFromB.size(); a++){
        event[tausFromB[a]].statusPos();
      }
      if (!pythia.moreDecays()) continue; // Go ahead and decay the undecayed taus
      
      // Now do all other taus 
      // Turn on all modes of decay
      pythia.readString("15:onMode = on");
      for(unsigned a =0; a < tausNotFromB.size(); a++){
        event[tausNotFromB[a]].statusPos(); // Mark the taus not from B as needing decaying
      }
      if (!pythia.moreDecays()) continue;
    } // end of if(tauToMuOnly)
    
   // Now do all your selection requirements to determine if we keep the event:
    
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

    // order mu pt vector
    // std::sort(muPtVec.begin(),muPtVec.end());
    // if (muPtVec[0] < 10 || muPtVec[1] < 10 ) continue; // Skip if top 2 pt muons have pt < 10

    // if it gets to here, then we're happy with the event
    iEvent++;
    for (unsigned a = 0; a < muPtVec.size(); a++){
      muPt.fill(muPtVec.at(a));
    }
    nMuInEvent.fill(nMuPos + nMuNeg);

    // Output the event to screen
    if (outputEvent)
      event.list();

    if (wanted && writeToHEPMC){
      // Construct new empty HepMC event and fill it.
      // Units will be as chosen for HepMC build, but can be changed
      // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)  
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      ToHepMC.fill_next_event( pythia, hepmcevt );

      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;
    }
  } //end of events loop

  // Statistics. Histograms. 
  pythia.stat();
  cout << muPt << nMuInEvent << endl;
  cout << "Number of events with pair: " << nWithPair << endl;

  // Done. 
  return 0;
}
