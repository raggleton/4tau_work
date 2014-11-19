#ifndef CUTS_H
#define CUTS_H

// For sorting methods, etc
#include "commonFunctions.h"

using std::cout;
using std::endl;

/**
 * Lots of fucntions that do cuts
 * Really this should be an object, but that requires more thinking...
 *
 * ALL DISTANCES IN MM, ALL MOMENTUM/ENERGY IN GeV
 *
 */

/**
 * This function takes two vectors of muons, and a muon test function,
 * and returns a std::pair (pT descending) if highest pT muons passing test function.
 * Used a lot in cutFlow program
 * @param  muons17toInf  std::vector of muons with pT > 17 GeV, descending pT order
 * @param  muons10to17   std::vector of muons with 17 > pT > 10 GeV, descending pT order
 * @param  muonCheck    function to test muons against
 * @return   std::pair of highest-pT muons passing cuts.
 */
std::pair<Track*, Track*> testMuons(std::vector<Track*> muons17toInf,
		  std::vector<Track*> muons10to17,
		  bool (*muonCheck)(Track*,Track*)) {

	Track *mu1(nullptr), *mu2(nullptr);
	bool foundMuonPair = false;
	std::vector<Track*>::iterator muA = muons17toInf.begin();
	while(!foundMuonPair && muA != muons17toInf.end()){

		// Need to make pairs among the 17toInf vector also, if size >= 2
		auto muB = muA;
		muB++;
		for (; muB != muA, muB != muons17toInf.end(); muB++) {
			if (*muA != *muB){
				if (muonCheck(*muA, *muB)) {
					mu1 = *muA;
					mu2 = *muB;
					foundMuonPair = true;
					break;
				}
			}
		}

		if(!foundMuonPair) {
			for (auto muB : muons10to17) {
				if (muonCheck(*muA, muB)) {
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

/**
 * This is like the above function but adds an extra argument for specifying
 * the deltaR(mu-mu) cut for the muonCheck function.
 * @param  muons17toInf  std::vector of muons with pT > 17 GeV, descending pT order
 * @param  muons10to17   std::vector of muons with 17 > pT > 10 GeV, descending pT order
 * @param  muonCheck    function to test muons against, which is dependent on deltaR param
 * @param  deltaR        deltaR(mu-mu) cut to use
 * @return   std::pair of highest-pT muons passing cuts.
 */
std::pair<Track*, Track*> testMuons(std::vector<Track*> muons17toInf,
		  std::vector<Track*> muons10to17,
		  bool (*muonCheck)(Track*, Track*, double),
		  double deltaR) {

	Track *mu1(nullptr), *mu2(nullptr);
	bool foundMuonPair = false;
	std::vector<Track*>::iterator muA = muons17toInf.begin();
	while(!foundMuonPair && muA != muons17toInf.end()){

		// Need to make pairs among the 17toInf vector also, if size >= 2
		auto muB = muA;
		muB++;
		for (; muB != muA, muB != muons17toInf.end(); muB++) {
			if (*muA != *muB){
				if (muonCheck(*muA, *muB, deltaR)) {
					mu1 = *muA;
					mu2 = *muB;
					foundMuonPair = true;
					break;
				}
			}
		}

		if(!foundMuonPair) {
			for (auto muB : muons10to17) {
				if (muonCheck(*muA, muB, deltaR)) {
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

/**
 * Calculates the transverse impact parameter, d_xy, as Delphes doesn't do it
 * unless you use the ImpactParameterSmearing module
 * @param  cand Track object
 * @return      Transverse impact parameter
 */
float calcDxy(Track* cand) {


	// calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    float xd =  cand->Xd;
    float yd =  cand->Yd;

    TLorentzVector candidateMomentum = cand->P4();

    // calculate impact paramater (after-smearing)
    float ang_mom = (xd*candidateMomentum.Py() - yd*candidateMomentum.Px());
    float dxy = ang_mom/cand->PT;

    cand->Dxy = dxy;
    return dxy;
}

/////////////////////////
// LOTS OF 2-MUON CUTS //
/////////////////////////

// Template these?

bool checkMuonsEta(Track* muA, Track* muB){
	if ((fabs(muA->Eta) < 2.1) && (fabs(muB->Eta) < 2.1)){
		return true;
	} else {
		return false;
	}
}

/**
 * Check to see if muA and muB satisfy muon IP criteria
 * @param  muA [description]
 * @param  muB [description]
 * @return     [description]
 */
bool checkMuonsIP(Track* muA, Track* muB){
	if( (fabs(muA->Zd) < 1.) // dZ < 0.1cm
		&& (fabs(muB->Zd) < 1.) // dZ < 0.1cm
		&& (fabs(calcDxy(muA)) < 0.3 ) // d0 < 0.03cm
		&& (fabs(calcDxy(muB)) < 0.3 ) // d0 < 0.03cm
	) {
		return true;
	} else {
		return false;
	}
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
	if (muA->PT > 17. && muB->PT > 10.
		&& checkMuonsIP(muA, muB)
		// && (fabs(muA->Zd) < 1.) // dZ < 0.1cm
		// && (fabs(muB->Zd) < 1.) // dZ < 0.1cm
		// // && (fabs(muA->Dxy) < 0.3 ) // d0 < 0.03cm
		// && (fabs(calcDxy(muA)) < 0.3 ) // d0 < 0.03cm
		// // && (fabs(muB->Dxy) < 0.3 ) // d0 < 0.03cm
		// && (fabs(calcDxy(muB)) < 0.3 ) // d0 < 0.03cm
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
 * These checks to see if muA and muB satisfy muon SS condtions
 * for signal region
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsSS(Track* muA, Track* muB){
	if (muA->Charge == muB->Charge){
		return true;
	} else {
		return false;
	}
}

/**
 * Check to see if muons pass IP and SS cuts
 * @param  muA [description]
 * @param  muB [description]
 * @return     [description]
 */
bool checkMuonsIPSSEta(Track* muA, Track* muB) {
	return checkMuonsEta(muA, muB) && checkMuonsIP(muA, muB) && checkMuonsSS(muA, muB);
}

/**
 * Check to see if muons pass IP and SS and deltaR cuts
 * @param  muA [description]
 * @param  muB [description]
 * @param  deltaR [description]
 * @return     [description]
 */
bool checkMuonsIPSSDREta(Track* muA, Track* muB, double deltaR) {
	if (checkMuonsEta(muA, muB)
		&& checkMuonsIP(muA, muB)
		&& checkMuonsSS(muA, muB)
		&& (muA->P4().DeltaR(muB->P4()) > deltaR)) {
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
	if (checkMuonsPTSS(muA, muB) && checkMuonsEta(muA, muB)){
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


/**
 * These checks to see if muA and muB satisfy all muon condtions,
 * but deltaR(mu-mu) > 1, used for sideband.
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsDR1(Track* muA, Track* muB){
	return checkMuons(muA, muB, 1.0);
}


/**
 * These checks to see if muA and muB satisfy all muon condtions
 * for signal region (dR(mu-mu) > 2)
 * @param  muA Higher pT muon
 * @param  muB Lesser pT muon
 * @return    TRUE if muA and muB pass cuts, FALSE otherwise
 */
bool checkMuonsSignal(Track* muA, Track* muB){
	return checkMuons(muA, muB, 2.0);
}


////////////////
// TRACK CUTS //
////////////////

/**
 * Check to see if track passes loose pT requirement (> 1 GeV)
 * @param  candTk Pointer to track object
 * @return        TRUE if track passes cuts, FALSE otherwise
 */
bool checkTrackPTLoose(Track* candTk){
	if (candTk->PT > 1.){
		return true;
	} else {
		return false;
	}
}

/**
 * Check to see if track passes loose IP requirement (dz < 1 cm, d0 < 1 cm)
 * (Delphes does all distacnes in mm)
 * @param  candTk Pointer to track object
 * @return        TRUE if track passes cuts, FALSE otherwise
 */
bool checkTrackIPLoose(Track* candTk){
	if ((fabs(candTk->Zd) < 10.) // dz < 1cm
		&& (fabs(calcDxy(candTk)) < 10.)){ // d0 impact parameter < 1cm
		return true;
	} else {
		return false;
	}
}

/**
 * Check to see if track passes Tight pT requirement (> 2.5 GeV)
 * (Delphes does all distacnes in mm)
 * @param  candTk Pointer to track object
 * @return        TRUE if track passes cuts, FALSE otherwise
 */
bool checkTrackPTTight(Track* candTk){
	if (candTk->PT > 2.5){
		return true;
	} else {
		return false;
	}
}

/**
 * Check to see if track passes tight IP requirement (dz < 0.04 cm, d0 < 0.02cm)
 * (Delphes does all distacnes in mm)
 * @param  candTk Pointer to track object
 * @return        TRUE if track passes cuts, FALSE otherwise
 */
bool checkTrackIPTight(Track* candTk){
	if ((fabs(candTk->Zd) < 0.4) // dz < 0.04cm
		&& (fabs(calcDxy(candTk)) < 0.2)) { // d0 < 0.02cm
		return true;
	} else {
		return false;
	}
}

/**
 * Check to see if track passes eta cut (|eta| < 2.4)
 * @param  candTk Pointer to track object
 * @return        TRUE if track passes cuts, FALSE otherwise
 */
bool checkTrackEta(Track* candTk) {
	if (fabs(candTk->Eta) < 2.4) {
		return true;
	} else {
		return false;
	}
}

/**
 * Check track against loose pT cuts, loose IP cut. No eta cut
 * - that's done in fillTrackVector since we have to rescale the track eta.
 * @param  candTk Pointer to track object
 * @return        TRUE if track passes cuts, FALSE otherwise
 */
bool checkTrackLoose(Track* candTk) {
	if (checkTrackPTLoose(candTk)
		&& checkTrackIPLoose(candTk)
		&& checkTrackEta(candTk)
		) {
		return true;
	} else {
		return false;
	}
}

/**
 * Check track against tight pT cuts, Tight IP cut, and eta cut/
 * @param  candTk Pointer to track object
 * @return        TRUE if track passes cuts, FALSE otherwise
 */
bool checkTrackTight(Track* candTk) {
	if (checkTrackPTTight(candTk)
		&& checkTrackIPTight(candTk)
		&& checkTrackEta(candTk)) {
		return true;
	} else {
		return false;
	}
}


/////////////////////
// TRACK-MUON CUTS //
/////////////////////

/**
 * Calculates deltaR(tk-mu) between candTk and mu1, and candTk and mu2.
 * If either is < deltaR param (default = 0.5) then fill the relevant vector
 * with pointer to that track.
 * @param candTk Track to compare against muons
 * @param mu1    One muon
 * @param mu2    Other muon
 * @param tks1   vector of track pointers, where dR(tk-mu1) < deltaR
 * @param tks2   vector of track pointers, where dR(tk-mu2) < deltaR
 * @param deltaR Optional deltaR cut for tk-mu. Default = 0.5
 */
void fillTrackVectors(Track* candTk, Track* mu1, Track* mu2,
				 std::vector<Track*>* tks1, std::vector<Track*>* tks2,
				 double deltaR_mutk = 0.5) {

	double dR1 = (candTk->P4()).DeltaR(mu1->P4());
	double dR2 = (candTk->P4()).DeltaR(mu2->P4());

	// do cut on track eta
	if (checkTrackEta(candTk)) {
		if (dR1 < deltaR_mutk) {
			tks1->push_back(candTk);
		}
		if (dR2 < deltaR_mutk) {
			tks2->push_back(candTk);
		}
	}
}

/**
 * @brief Rescales track eta & phi using muon eta & phi
 * @details  DON'T set using the TLorentzVEctor from P4() method
 * - it won't store the new Eta/Phi
 * Instead do candTk->Eta = ... etc
 * Then if you call candTk->P4() again, it'll give the new updated TLorentzVector
 * Also, use the DelptaPhi function to calculate the difference in phi correctly,
 * otherwise you'll get lost in factors of 2pi
 * @param tk Track to rescale
 * @param mu Muon to use for rescaling
 * @param alpha Rescale factor. 1 = off (no rescaling)
 */
void rescaleTk(Track* tk, Track* mu, double rescaleFactor) {
	if (rescaleFactor != 1) { // to save time?
		tk->Eta = mu->Eta + (rescaleFactor * (tk->Eta - mu->Eta));
		tk->Phi = mu->Phi + (rescaleFactor * (tk->P4()).DeltaPhi(mu->P4()));
	}
}

/**
 * Calculates deltaR(tk-mu) between candTk and mu1, and candTk and mu2.
 * If either is < deltaR param (default = 0.5) then fill the relevant vector
 * with pointer to that track.
 * Also does rescaling of track eta and phi, and checks track eta is in barrel
 * @param candTk Track to compare against muons
 * @param mu1    One muon
 * @param mu2    Other muon
 * @param tks1   vector of track pointers, where dR(tk-mu1) < deltaR
 * @param tks2   vector of track pointers, where dR(tk-mu2) < deltaR
 * @param rescale rescale value (alpha) for scaling MC to data. 1 = no scaling
 * @param deltaR Optional deltaR cut for tk-mu. Default = 0.5
 */
void fillTrackVectorsWithRescale(Track* candTk, Track* mu1, Track* mu2,
				 std::vector<Track*>* tks1, std::vector<Track*>* tks2,
				 double rescaleFactor, double deltaR_mutk = 0.5) {

	// store pre-scaled values
	double dR1 = (candTk->P4()).DeltaR(mu1->P4());
	double dR2 = (candTk->P4()).DeltaR(mu2->P4());

	// Skip if nether track will pass dR selection for either muon
	if (((dR1*rescaleFactor) > deltaR_mutk)
		&& ((dR2*rescaleFactor) > deltaR_mutk)) {
		return;
	}

	if (dR1 < dR2) {
		rescaleTk(candTk, mu1, rescaleFactor);
		if (checkTrackEta(candTk) && ((dR1*rescaleFactor) < deltaR_mutk)) {
			// cout << "prescaling: dR1: " << dR1 << ", expect " << rescaleFactor*dR1;
			// cout << " postscaling: dR1: " << (candTk->P4()).DeltaR(mu1->P4()) << " factor: " << (candTk->P4()).DeltaR(mu1->P4())/dR1 << endl;
			tks1->push_back(candTk);
		}
	} else {
		rescaleTk(candTk, mu2, rescaleFactor);
		if (checkTrackEta(candTk) && ((dR2*rescaleFactor) < deltaR_mutk)) {
			// cout << "prescaling: dR2: " << dR2 << ", expect " << rescaleFactor*dR2;
			// cout << " postscaling: dR2: " << (candTk->P4()).DeltaR(mu2->P4()) << " factor: " << (candTk->P4()).DeltaR(mu2->P4())/dR2 << endl;
			tks2->push_back(candTk);
		}
	}

}


/**
 * Check to see if track and muon have opposite charges
 * @param  candTk Track
 * @param  mu     Muon
 * @return        TRUE if opposite charge, FALSE otherwise
 */
bool checkTkMuOS(Track* candTk, Track* mu) {
	return ((candTk->Charge) * (mu->Charge)) < 0;
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