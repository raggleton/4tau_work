#ifndef CUTS_H
#define CUTS_H

// For sorting methods, etc
#include "commonFunctions.h"

using std::cout;
using std::endl;

/**
 * Lots of fucntions that do cuts
 * Really this should be an object, but that requires more thinking...
 */

/**
 * This function takes two vectors of muons, and a muon test function,
 * and returns a std::pair (pT descending) if highest pT muons passing test function.
 * Used a lot in cutFlow program
 * @param  muons17toInf  std::vector of muons with pT > 17 GeV, descending pT order
 * @param  muons10to17   std::vector of muons with 17 > pT > 10 GeV, descending pT order
 * @param  checkMuons    function to test muons against
 * @return   std::pair of highest-pT muons passing cuts.
 */
std::pair<Track*, Track*> testMuons(std::vector<Track*> muons17toInf, 
		  std::vector<Track*> muons10to17, 
		  bool (*checkMuons)(Track*,Track*)) {

	Track *mu1(nullptr), *mu2(nullptr);
	bool foundMuonPair = false;
	std::vector<Track*>::iterator muA = muons17toInf.begin();
	while(!foundMuonPair && muA != muons17toInf.end()){
		
		// Need to make pairs among the 17toInf vector also, if size >= 2
		auto muB = muA;
		muB++;
		for (; muB != muA, muB != muons17toInf.end(); muB++) {
			if (*muA != *muB){
				if (checkMuons(*muA, *muB)) {
					mu1 = *muA;
					mu2 = *muB;
					foundMuonPair = true;
					break;
				}
			}
		}

		if(!foundMuonPair) {
			for (auto muB : muons10to17) {
				if (checkMuons(*muA, muB)) {
					mu1 = *muA;
					mu2 = muB;
					foundMuonPair = true;
					break;
				}
			}
		}
		muA++;
	}
	std::pair<Track*, Track*> p(mu1, mu2);
	return p;
}

/////////////////////////
// LOTS OF 2-MUON CUTS //
/////////////////////////

// Template these?

/**
 * These checks to see if muA and muB satisfy all muon condtions,
 * but deltaR(mu-mu) > 1, used for sideband.
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsAllDR1(Track* muA, Track* muB){
	return checkMuons(muA, muB, 1.0);
}

/**
 * These checks to see if muA and muB satisfy all muon condtions 
 * for signal region (dR(mu-mu) > 2)
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsAllSignal(Track* muA, Track* muB){
	return checkMuons(muA, muB, 2.0);
}

/**
 * These checks to see if muA and muB satisfy muon pT condtions
 * & IP conditions
 * NB IP cuts only work with Track objects starting in Delphes 3.1.2
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsPT(Track* muA, Track* muB){
	if (muA->PT > 17. && muB->PT > 10
		&& (fabs(muA->Zd) < 1.) // dZ < 0.1cm
		&& (fabs(muB->Zd) < 1.) // dZ < 0.1cm
		&& (fabs(muA->Dxy) < 0.3 ) // d0 < 0.03cm
		&& (fabs(muB->Dxy) < 0.3 ) // d0 < 0.03cm
		){
		return true;
	} else {
		return false;
	}
}

/**
 * These checks to see if muA and muB satisfy muon pT and SS condtions 
 * for signal region
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsPTSS(Track* muA, Track* muB){
	if (checkMuonsPT(muA, muB) && (muA->Charge == muB->Charge)){
		return true;
	} else {
		return false;
	}
}

/**
 * These checks to see if muA and muB satisfy muon pT, SS and eta condtions 
 * for signal region
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsPTSSEta(Track* muA, Track* muB){
	if (checkMuonsPTSS(muA, muB) 
		&& (fabs(muA->Eta) < 2.1) && (fabs(muB->Eta) < 2.4)){
		return true;
	} else {
		return false;
	}
}

/**
 * These checks to see if muA and muB satisfy all muon condtions, 
 * with user-defined deltaR cut
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @param  deltaR dR(mu-mu) cut
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuons(Track* muA, Track* muB, double deltaR){
	if (checkMuonsPTSSEta(muA, muB)	&& (muA->P4().DeltaR(muB->P4()) > deltaR)){
		return true;
	} else {
		return false;
	}
}


/////////////////
// OLD VERSION //
/////////////////

/**
 * These checks to see if muA and muB satisfy all muon condtions
 * apart from impact params (as vars not stored)
 * (Old version kept incase you use RawGenMuons branch instead of GenMuons)
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @param  deltaR dR(mu-mu) cut
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuons(GenParticle* muA, GenParticle* muB, double deltaR){
	if (muA->PT > 17. && muB->PT > 10.
		&& (muA->Charge == muB->Charge)
		&& (fabs(muA->Eta) < 2.1)
		&& (fabs(muB->Eta) < 2.1)
		&& ((muA->P4().DeltaR(muB->P4())) > deltaR)
		){
		return true;
	} else {
		return false;
	}
}

#endif