#ifndef MYHOOKS_H
#define MYHOOKS_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

bool checkBorC(int ID) {
  if (abs(ID) == 4 || abs(ID) == 5)
    return true;
  else
    return false;
}

// Write own derived UserHooks class.
// This User Hook vetoes quark-gluon scatter events 
// unless the incoming quark is c or b,
// and vetoes unless gluon decays to ccbar or bbbar
class QGScatterHook : public UserHooks {

  public:

    // Constructor 
    QGScatterHook(bool verbose = true): 
      DEBUG(verbose)
    {
      // any ctor stuff here
    }

    // Destructor 
    ~QGScatterHook() {
      // any dtor stuff here
    }

    // Allow a veto for the hard process (and associated resonance decays) 
    virtual bool canVetoProcessLevel() { return true; }

    // At this stage, the process record typically contains the two beams in 
    // slots 1 and 2, the two incoming partons to the hard process in slots 
    // 3 and 4, the N (usually 1, 2 or 3) primary produced particles in slots 5 
    // through 4 + N, and thereafter recursively the resonance decay chains, 
    // if any.
    virtual bool doVetoProcessLevel(Event& process) {
      
      // 1. Want qg scatter to only allow cg or bg incoming
      // if (process[3].id() == 21) {
      //   if (process[4].idAbs() == 4 || process[4].idAbs() == 5) {
      //     if (DEBUG) {
      //       cout << "doVetoProcessLevel PDGIDs" << endl;
      //       cout << "slot 1: " << process[1].id() << ", slot 2: " << process[2].id() 
      //       << ", slot 3: " << process[3].id() << ", slot 4: " << process[4].id() 
      //       << ", slot 5: " << process[5].id() << ", slot 6: " << process[6].id() 
      //       << ", slot 7: " << process[7].id() << endl;
      //       cout << "Gluon daughters index: " << process[4].daughter1() << " to " << process[4].daughter2() << endl;
      //     }
      //     return false;
      //   } else {
      //     return true;
      //   }
      // } else {
      //   if (process[3].idAbs() == 4 || process[3].idAbs() == 5) {
      //     if (DEBUG) {
      //       cout << "doVetoProcessLevel PDGIDs" << endl;
      //       cout << "slot 1: " << process[1].id() << ", slot 2: " << process[2].id() 
      //       << ", slot 3: " << process[3].id() << ", slot 4: " << process[4].id() 
      //       << ", slot 5: " << process[5].id() << ", slot 6: " << process[6].id() 
      //       << ", slot 7: " << process[7].id() << endl;
      //       cout << "Gluon daughters index: " << process[3].daughter1() << " to " << process[3].daughter2() << endl;
      //     }
      //     return false;
      //   } else {
      //     return true;
      //   }
      // }

      // 2. Want the gluon to only decay to ccbar or bbbar
      // Since the continuing quark doesn't change flav, 
      // we require all 3 particles from process[5] to process[7] must be b/c
      if (checkBorC(process[5].idAbs()) && checkBorC(process[6].idAbs()) && checkBorC(process[7].idAbs())){
        if (DEBUG) { cout << "slot 1: " << process[1].id() << ", slot 2: " << process[2].id() 
            << ", slot 3: " << process[3].id() << ", slot 4: " << process[4].id() 
            << ", slot 5: " << process[5].id() << ", slot 6: " << process[6].id() 
            << ", slot 7: " << process[7].id() << endl; }
        return false;
      } else {
        return true;
      }
    }   

  private:
    bool DEBUG;
};

#endif