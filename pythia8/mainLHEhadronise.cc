// main12.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This takes in a LHE file and hadronises it, and output to hepmc
// Based on main12.cc in pythia/examples

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"   
#include "HepMC/IO_GenEvent.h"

using namespace Pythia8; 

/**
 * This method analyses the decay products from a tau. Has to go through the complete decay chain
 * until it finds all the final state particles. Then it figure sout which are charged, and which are muons.
 * nProngs: returns the number of charged particles from the tau (includes muons!)
 * nMu: returns the number of stable muons from the tau
 * current: the immediate decays products of the tau, passed to event from method that calls this fn.
 */
void lookAtTauProducts(Event& event, int &nProngs, int &nMu, std::vector<int> current) {
  std::vector<int> history; // hold all unique particles in the decay chain (stores event posiiton number)
  // std::vector<int> current; // holds position no.s for current step
  std::vector<int> next; // holds decay products, not nec. all unique

  while (current.size()>0){ // if current > 0 we haven't exhausted all the particles
    for (unsigned a = 0; a < current.size(); a++){
      // Check 1 - is this already in current?
      // Could probably do more efficiently using the unique function on std::vector
      bool alreadyDone = false;

      for (unsigned b = 0; b < a; b++){
        if ((current[a] == current[b]) && (a != 0)) {
          alreadyDone = true;
        }
      }

      // Check 2 - is this already in history?
      if (!alreadyDone){
        for (unsigned c = 0; c < history.size(); c++){
          if ((current[a] == history[c]) && (c!=0)){
            alreadyDone = true;
          }
        }
      }

      // Check 3 - is this final state?
      if (!alreadyDone){
        if (event[current[a]].status() > 0) {
          // Check if muon
          if (event[current[a]].idAbs() == 13){
            nMu++;
            nProngs++;
          } else {
            // Check if charged
            if (event[current[a]].isCharged()){
              nProngs++;
            }
          }
          // Load it into history
          history.push_back(current[a]);
        } else {
          // Load its daughters no. into next
          if (event[current[a]].daughter2() == 0) // if only 1 daughter, the second daughter will return 0
            next.push_back(event[current[a]].daughter1());
          for (int d = event[current[a]].daughter1(); d <= event[current[a]].daughter2(); d++) {
            next.push_back(d);
          } 
          
          // Load it into history
          history.push_back(current[a]);
        }
      }// end of alreadyDone
    } // end of loop over current

    // Clear current - don't need to do as = assignment auto does this
    // current.clear();
    // Copy next into current
    current = next;
    // Empty next
    next.clear();
  } // end of while(!done)

}


int main(int argc, char* argv[]) {

  // bool outputEvent  = true; // output entire event listing to STDOUT (long!), for debugging only
  bool writeToHEPMC = true; // output to HEPMC

  // Check that correct number of command-line arguments
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments. \n "
         << " You are expected to provide one output file name. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // cout << "Outputting to " << argv[1] << endl;

  // Interface for conversion from Pythia8::Event to HepMC event. 
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[1], std::ios::out);

  //  Number of listed events. Allow for possibility of a few faulty events.
  int nPrintLHA  = 1;             
  int nPrintRest = 0;             
  int nAbort     = 10;
  int nMaxEvent  = 500000;
  
  // Generator           
  Pythia pythia;                            
  Event& event = pythia.event;

  // Do random seed
  //  a value 0 gives a random seed based on the time
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");

  // No automatic event listings - do it manually below.
  pythia.readString("Next:numberShowLHA = 0"); 
  pythia.readString("Next:numberShowInfo = 0"); 
  pythia.readString("Next:numberShowProcess = 0"); 
  pythia.readString("Next:numberShowEvent = 0"); 

  // Initialize Les Houches Event File run.
  pythia.readString("Beams:frameType = 4"); // the beam and event information is stored in a Les Houches Event File
  pythia.readString("Beams:LHEF = ../Signal_1prong_500K_bare/GG_H_aa_8_4taus_decay_500K_2-single.lhe");   
  // pythia.readString("Beams:LHEF = reduced_GG_H_aa_4taus_2.lhe");   
  
  // pythia.readString("ProcessLevel:all = off");   
  // pythia.readString("PartonLevel:all = off");   
  // pythia.readString("HadronLevel:all = off");   
  pythia.init();   

  // Some basic historgrams
  Hist nMuInEvent("number of muons in an event", 10, -0.5, 9.5); 
  Hist muPt("pT muons in an event", 40, 0.0, 20.0); 
  
  // Set counters.
  int iPrintLHA  = 0;             
  int iPrintRest = 0;             
  int iAbort     = 0;
  int iFile      = 1;

  int nWanted(0);

  // Begin event loop   
  for (int iEvent = 0; ; ++iEvent) {
    bool wanted = false;

    // Generate until none left in input file or get to nMaxEvent
    if (!pythia.next() || iEvent > nMaxEvent) {
      if (++iAbort < nAbort) continue;
      break;
    }

    // cout << "****Event" << endl;
  
    // List first few Les Houches and other events.
    if (pythia.info.isLHA() && iPrintLHA < nPrintLHA) {     
      pythia.LHAeventList();               
      pythia.info.list();          
      pythia.process.list();          
      pythia.event.list();  
      ++iPrintLHA;         
    } else if (!pythia.info.isLHA() && iPrintRest < nPrintRest) {     
      pythia.info.list();          
      pythia.process.list();          
      pythia.event.list();           
      ++iPrintRest;         
    }                 

    // Look for 1 prong taus only
    int n1Prong(0), n3Prong(0), nMus(0); //xProng includes muons
    // Look for muons among decay products (also from charm/tau/...).
    std::vector<double> muPtVec;

    // Loop over all particles in event
    for (int i = 0; (i < event.size()) && (n3Prong==0); ++i) {
      if ((event[i].idAbs() == 15 ) && (event[i].status() == -22)){ // loop over all taus from a_0

        int nProngs(0), nMu(0);
        std::vector<int> daughters;
        for (int d = event[i].daughter1(); d <= event[i].daughter2() && d>0; d++){
          daughters.push_back(d);
        }
        lookAtTauProducts(event, nProngs, nMu, daughters);
        // cout << "There were " << nMu << " muons from this tau, and " << nProngs << " charged tracks from this tau" << endl;
        if (nProngs == 1) n1Prong++;
        else if (nProngs > 1) n3Prong++;
        nMus += nMu;
      }


      int id = event[i].id();  
      if ((abs(id) == 13) && (event[i].status() > 0)){ 
          muPtVec.push_back(event[i].pT());
      }
      

    } // end loop over particles in event
    
    // After looking at all taus, have we got 2 muons and 2 1-prong decays? 
    // (the muPtVec.size() test is kinda redundant, 
    // as the 2 gen muons should provide 2 muons in that vec)
    if ((n1Prong == 4) && (nMus>=2) && (muPtVec.size()>1)){

      // order mu pt vector
      std::sort(muPtVec.begin(),muPtVec.end(), std::greater<int>());

      // Emulate HLT - HLT_Mu17_Mu8
      if (muPtVec[0]>17 && muPtVec[1]>8 ) {
        wanted = true;
        nWanted++;
        for (unsigned a = 0; a < muPtVec.size(); a++){
          muPt.fill(muPtVec.at(a));
        }
        nMuInEvent.fill(muPtVec.size());

      // cout << "+++++++++++++++++ YAYYYY +++++++++++++" <<endl;
      }
    } else {
      wanted = false;
    }

    // cout << "This event had " << n1Prong << " 1-prong taus, " << n3Prong << " 3-prong taus and " << nMus << " tau to mu decays." << endl;
    if (n3Prong+n1Prong > 4) {
      cout << "*********** SHIT > 4 --------------------------------------" <<endl;
            // pythia.event.list();  
    }

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

    muPtVec.clear();
  
  } // End of event loop.        
  cout << " Number of useful events: " << nWanted << "/" << nMaxEvent << endl;
  // Give statistics. Print histogram.
  pythia.stat();
  cout << muPt << nMuInEvent << endl;

  // Done.                           
  return 0;
}
