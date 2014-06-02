#ifndef CUTS_H
#define CUTS_H

// For sorting methods, etc
#include "commonFunctions.h"

using std::cout;
using std::endl;


/**
 * Object that does cuts for all programs
 */
// class Cuts
// {
// 	private:
// 		std::vector<T*> muons;
// 	public:
// 		// Ctor
// 		Cuts():
// 		{

// 		}
// };

std::pair<GenParticle*, GenParticle*> testMuons(std::vector<GenParticle*> muons17toInf, 
		  std::vector<GenParticle*> muons10to17, 
		  bool (*checkMuons)(GenParticle*,GenParticle*)) {

	GenParticle *mu1(nullptr), *mu2(nullptr);
	bool foundMuonPair = false;
	// if (!p.first) cout << "p.first NULL" << endl;
	std::vector<GenParticle*>::iterator muA = muons17toInf.begin();
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
	std::pair<GenParticle*, GenParticle*> p(mu1, mu2);
	// if (!p.first) cout << "p.first NULL" << endl; else cout << "p.first not null" << endl;
	return p;
}

/**
 * These checks to see if muA and muB satisfy all muon condtions,
 * but deltaR(mu-mu) > 1, used for sideband.
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsAllDR1(GenParticle* muA, GenParticle* muB){
	if (muA->PT > 17. && muB->PT > 10
		&& (muA->Charge == muB->Charge)
		&& (fabs(muA->Eta) < 2.1)
		&& (fabs(muB->Eta) < 2.4)
		&& ((muA->P4().DeltaR(muB->P4())) > 1.)
		){
		return true;
	} else {
		return false;
	}
}

/**
 * These checks to see if muA and muB satisfy all muon condtions 
 * for signal region
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsAllSignal(GenParticle* muA, GenParticle* muB){
	if (muA->PT > 17. && muB->PT > 10
		&& (muA->Charge == muB->Charge)
		&& (fabs(muA->Eta) < 2.1)
		&& (fabs(muB->Eta) < 2.4)
		&& ((muA->P4().DeltaR(muB->P4())) > 2.)
		){
		return true;
	} else {
		return false;
	}
}

/**
 * These checks to see if muA and muB satisfy muon pT condtions
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsPT(GenParticle* muA, GenParticle* muB){
	if (muA->PT > 17. && muB->PT > 10){
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
bool checkMuonsPTSS(GenParticle* muA, GenParticle* muB){
	if (muA->PT > 17. && muB->PT > 10
		&& (muA->Charge == muB->Charge)
		){
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
bool checkMuonsPTSSEta(GenParticle* muA, GenParticle* muB){
	if (muA->PT > 17. && muB->PT > 10
		&& (muA->Charge == muB->Charge)
		&& (fabs(muA->Eta) < 2.1)
		&& (fabs(muB->Eta) < 2.4)
		){
		return true;
	} else {
		return false;
	}
}


#endif