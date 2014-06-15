// Based off main15.cc in examples folder
// 
// This creates QCD b samples, which can optionally decay to muons

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
  // A value 0 gives a random seed based on the time
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  
  ////////////////////////
  // Number of events.  //
  ////////////////////////
  // For HLT. NoHLT has about 60X HLT amount (320K NoHLT evt for 5K HLT evnt)
  // Warning, 5K events ~900MB hepmc file and takes ~5 min.
  // Warning, 50K events ~9GB hepmc file and takes ~40 min.
  int nEvent = 500;
  pythia.readString("Next:numberShowEvent = 00");
  // pythia.readString("Next:numberShowProcess = 100");
  

  ///////////////////////////
  // Setup for pp at 8 TeV //
  ///////////////////////////
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 8000");

  // Making it like PythiaUEZ2starSettings
  pythia.readString("PDF:pSet = 8"); // CTEQ6L1

  ////////////////////////////////////////////////////
  // Simulate b production above given pTmin scale. //
  ////////////////////////////////////////////////////
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

  ////////////////////////////////////////////////
  // Set B hadron decay modes to mu or tau only //
  ////////////////////////////////////////////////

  // All B hadrons from PDG
  // NOTE: vector only contains the particle IDs *NOT* the anti particle IDs ( = -PDGID )!!
  // Deal with this below
   std::vector<int> bCodes{ 511,521,10511,10521,513,523,10513,10523,20513,20523,515,525,531,10531,533,10533,
    20533,535,541,10541,543,10543,20543,545,51,10551,100551,110551,200551,210551,553,10553,20553,30553,
    100553,110553,120553,130553,200553,210553,220553,300553,10860,9000553,11020,9010553,555,10555,20555,
    100555,110555,120555,200555,557,100557,5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,5314,
    5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,5514,5524,5532,5534,
    5542,5544,5554};
  int nCodes = bCodes.size();

  if (muOnly){
    // Set B hadrons to decay to modes involving a muon or tau
    for (int iC = 0; iC < nCodes; ++iC) {
      // Check PDGID is in PYTHIA
      if(! pythia.particleData.isParticle(bCodes[iC])) continue;
      
      // Get particle name.
      // If excited state, then just skip it, we want it to decay to less excited states.
      if (pythia.particleData.name(bCodes[iC]).find("*") != std::string::npos) continue;
      if (pythia.particleData.name(-bCodes[iC]).find("*") != std::string::npos) continue;
      
      // Turn off all decay modes first (for particle & antiparticle)
      std::stringstream sstm;
      sstm << bCodes[iC] << ":onMode = off";
      std::string command = sstm.str();
      pythia.readString(command);
      pythia.readString("-"+command);
      sstm.str("");
      
      // Now just turn on tau or mu ones
      sstm << bCodes[iC] << ":onIfAny = 13 15";
      command = sstm.str();
      pythia.readString(command);
      pythia.readString("-"+command);
    }
    pythia.particleData.list(511);
    pythia.particleData.list(-511);
  }

  // For tau, turn off decays. We want any tau from B hadrons to decay to muons,
  // but all other taus can decay however they want.
  if(tauToMuOnly)
    pythia.readString("15:onMode = off");
    pythia.readString("-15:onMode = off");
  
  std::vector<int> tausFromB;
  std::vector<int> tausNotFromB;

  // Initialize 
  pythia.init();

    pythia.particleData.list(511);
    pythia.particleData.list(-511);


  // Some basic histograms
  Hist nMuInEvent("number of muons in an event (HLT)", 10, -0.5, 9.5); 
  Hist muPt("pT muons in an event (HLT)", 40, 0.0, 40.0); 
  Hist muPtNoHLT("pT muons in an event (NoHLT)", 40, 0.0, 40.0); 
  int nWithMuPair = 0;

  ///////////////////////
  // Begin event loop. //
  ///////////////////////
  int iEvent      = 0; // counts events passing HLT, NOT NoHLT
  int lastiEvent  = 0; // for outputting to screen every outputEvery events
  int outputEvery = 50; 

  while(iEvent < nEvent) {
    bool wantedHLT = false;
    bool wantedNoHLT = false;

    if ((iEvent % outputEvery == 0) && (iEvent!= lastiEvent)){
      lastiEvent = iEvent;
      cout << "iEvent: " << iEvent << endl;
    }

    if (tauToMuOnly){
      // have to turn off tau decays for each event
      pythia.readString("15:onMode = off"); 
      pythia.readString("-15:onMode = off"); 
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
      ++nWithMuPair;
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
  cout << "Number of events with pair, & passing HLT: " << nWithMuPair << endl;

  // Done. 
  return 0;
}
