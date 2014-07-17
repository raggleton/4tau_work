#include "commonFunctions.h"
#include "cuts.h"

using std::cout;
using std::endl;

/**
 * Outputs # of events passing each stage of cuts
 * Not very sophisticated, could be made better!
 */

// Yeah, I know global vars are hated, but it's handy for quickly debugging fns
bool DEBUG = false;

// For sorting track vectors by pT
void sortTrackVector(std::vector<Track*>& tk){
	std::sort(tk.begin(), tk.end(), sortByPT<Track>);
}

/**
 * This checks that the pair p has non-null pointers (i.e. 2 muons passed cuts)
 * If so, increment the iterator, and return TRUE
 * @param  p        std::pair of objects. Must be const as not changing the pair, only checking it, 
 *                  and allows us to do testResults(testMuons(...)) in one line.
 * @param  counter  reference to the counter we want to increment
 * @return    TRUE if pair is valid, FALSE otherwise
 */
bool testResults(const std::pair<Track*, Track*> &p, int &counter) {
	if (p.first && p.second) {
		if (DEBUG) cout << "p.first and p.second not null" << endl;
		counter++; 
		return true;
	} else { 
		return false;
	}
}

void cutFlow(int argc, char* argv[])
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	ProgramOpts pOpts(argc, argv);
	
	// MCsource source     = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	// bool doSignal       = pOpts.getSignal(); // for signal or not
	// bool doMu           = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = pOpts.getMuOrdering(); // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	// bool doHLT          = pOpts.getHLT(); // whether to use MC that has HLT cuts already applied or not.
	DEBUG = pOpts.getVerbose();

	// Create chain of root trees
	TChain chain("Delphes");
	addInputFiles(&chain, &pOpts);
	pOpts.printProgramOptions();

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

	// Get pointers to branches used in this analysis
	// Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
	// and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
	// TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchTracks   = treeReader->UseBranch("Track");
	// TClonesArray *branchGenMuons = treeReader->UseBranch("OnlyGenMuons");
	TClonesArray *branchGenMuons = treeReader->UseBranch("GenMuon");
	// TClonesArray *branchStable   = treeReader->UseBranch("StableParticle");
	// TClonesArray *branchAll      = treeReader->UseBranch("AllParticle");

	//////////////////////
	// Loop over events //
	//////////////////////
	Long64_t numberOfEntries = getNumberEvents(treeReader, &pOpts);
	cout << "Running over " << numberOfEntries << " events" << endl;

	bool stop = false; // used to stop the loop, for debugging/testing

	std::vector<int> cutCount(8); // hold # evts passing cut
	std::vector<int>::iterator it = cutCount.begin();

	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		if (DEBUG) cout << "*** Event" << endl;
		it = cutCount.begin();

		if (DEBUG) cout << "Testing if >=2 muons" << endl;
		if (branchGenMuons->GetEntries() >= 2) {
			(*it)++; 
		} else continue; // skip if <2 muons!
		it++;

		{
			//////////////////////////////////////////////////////
			// Get the hard interaction particles for signal MC //
			// No selection cuts applied (only >=2 muons)       //
			//////////////////////////////////////////////////////
			
			// if (doSignal) {
			// 	GenParticle *a1(nullptr), *a2(nullptr);
			// 	// Get a0s
			// 	for(int j = 0; j < branchAll->GetEntries(); j++){
			// 		GenParticle *candHa = (GenParticle*) branchAll->At(j);
			// 		// cout << j << " ID: " << candHa->PID << " status: " << candHa->Status << endl;
				
			// 		// if (fabs(candHa->PID)== 12|| fabs(candHa->PID)== 14|| fabs(candHa->PID)==16 ){
			// 		// 	histNuPt->Fill(candHa->PT);
			// 		// } 

			// 		// Which is "1" and "2" is arbitrary here.
			// 		if ((fabs(candHa->PID)==36) && (fabs(candHa->Status)==62)){
			// 			if (a1==0){
			// 				a1=candHa;
			// 			} else {
			// 				a2=candHa;
			// 			}
			// 		}
			// 	}

			// 	// Get the tau daughters from a1 and a2 (no pT ordering)
			// 	GenParticle *tau1a(nullptr), *tau1b(nullptr), *tau2a(nullptr), *tau2b(nullptr);
			// 	tau1a = (GenParticle*) branchAll->At(a1->D1);
			// 	tau1b = (GenParticle*) branchAll->At(a1->D2);
			// 	tau2a = (GenParticle*) branchAll->At(a2->D1);
			// 	tau2b = (GenParticle*) branchAll->At(a2->D2);

			// 	// TLorentzVector tau1aMom,tau1bMom, tau2aMom, tau2bMom;
			// 	// tau1aMom = tau1a->P4();
			// 	// tau1bMom = tau1b->P4();
			// 	// tau2aMom = tau2a->P4();
			// 	// tau2bMom = tau2b->P4();

			// 	// histDRa1->Fill(tau1aMom.DeltaR(tau1bMom));
			// 	// histDRa2->Fill(tau2aMom.DeltaR(tau2bMom));
				
			// 	GenParticle *charged1a = getChargedObject(branchAll, tau1a);
			// 	GenParticle *charged1b = getChargedObject(branchAll, tau1b);
			// 	GenParticle *charged2a = getChargedObject(branchAll, tau2a);
			// 	GenParticle *charged2b = getChargedObject(branchAll, tau2b);
				
			// 	// This selects events where each tau only has 1 charged product...
			// 	// dunno what to do about evts where the tau decays into charged things *including* muon
			// 	if (charged1a && charged1b && charged2a && charged2b){
					
			// 		// To hold mu and tracks from each tau
			// 		GenParticle* muTruth1;
			// 		GenParticle* trackTruth1;
			// 		GenParticle* muTruth2;
			// 		GenParticle* trackTruth2;

			// 		// Assign charged products to be mu or track
			// 		bool truth1HasMu = assignMuonAndTrack(muTruth1, trackTruth1, *charged1a, *charged1b);				
			// 		bool truth2HasMu = assignMuonAndTrack(muTruth2, trackTruth2, *charged2a, *charged2b);

			// 		// NOTE: muons are NOT pT ordered

			// 		if (!truth1HasMu || !truth2HasMu) {
			// 			cout << "Problem, no truth mu for 1 and/or 2!" << endl;
			// 		} else { 
						
			// 			// Do m1 distribution in bins of m2 - for MC truth (is it actually correlated?)
			// 			double m1(0.);
			// 			double m2(0.);
						
			// 			// Assign m1 to higher pT muon
			// 			if (muTruth1->PT > muTruth2->PT) {
			// 				m1 = (muTruth1->P4()+trackTruth1->P4()).M();
			// 				m2 = (muTruth2->P4()+trackTruth2->P4()).M();
			// 			} else {
			// 				m2 = (muTruth1->P4()+trackTruth1->P4()).M();
			// 				m1 = (muTruth2->P4()+trackTruth2->P4()).M();
			// 			}

			// 			// Randomly swap trk-mu pairs 1<->2 if desired
			// 			if(swapMuRandomly){
			// 				double randNum = (double)rand() / RAND_MAX;
			// 				if (randNum > 0.5){
			// 					double tmp = m2;
			// 					m2 = m1;
			// 					m1 = tmp;
			// 				}
			// 			}

			// 			cout << m1 << "     " << m2 << endl;
			// 			// if(m2 < 1.)
			// 			// 	histM1_truth_0to1->Fill(m1);
			// 			// else if (m2 < 2.)
			// 			// 	histM1_truth_1to2->Fill(m1);
			// 			// else if (m2 < 3.)
			// 			// 	histM1_truth_2to3->Fill(m1);
			// 			// else
			// 			// 	histM1_truth_3toInf->Fill(m1);
			// 		}
			// 	} 
			// } // end if(doSignal)// MC truth stuff - not used
		}

		//////////////////////////////////////////////////////////////////////
		// Now, do more general MC stuff:                                   //
		// get the two highest pT muons in the event, store their pT        //
		// and pointers to the Tracks                                 //
		//////////////////////////////////////////////////////////////////////
		
		// Fill vectors with muons, based on pT
		std::vector<Track*> muons10to17;
		std::vector<Track*> muons17toInf;
		for (int i = 0; i < branchGenMuons->GetEntries(); i++){
			Track* cand = (Track*) branchGenMuons->At(i);
			if (cand->PT > 17) {
				muons17toInf.push_back(cand);
			} else if (cand->PT > 10) {
				muons10to17.push_back(cand);
			}
		}

		// Check to see if we can skip the event if not enough muons
		// Don't enable this as all HLT events so unnecessary
		// if (!(muons17toInf.size() >= 1 && (muons17toInf.size() + muons10to17.size()) >= 2)) continue;

		// Sort both vectors by descending pT
		std::sort(muons17toInf.begin(), muons17toInf.end(), sortByPT<Track>);
		std::sort(muons10to17.begin(), muons10to17.end(), sortByPT<Track>);

		// // Now randomly swap mu1 - mu2 if desired
		// Track *origMu1(nullptr), *origMu2(nullptr);
		// origMu1 = mu1;
		// origMu2 = mu2;
		// if (swapMuRandomly){
		// 	double randNum = (double)rand() / RAND_MAX;
		// 	// histRand->Fill(randNum);
		// 	if (randNum > 0.5){
		// 		mu1 = origMu2;
		// 		mu2 = origMu1;
		// 	}
		// }

		////////////////////
		// Muon selection //
		////////////////////
		if (DEBUG) cout << "Testing if mu with pT > 17" << endl;
		if (muons17toInf.size() > 0) {
			(*it)++;
		} else {
			continue;
		}
		it++;

		if (DEBUG) cout << "Testing if mu with pT > 10" << endl;
		if (muons10to17.size() > 0 || muons17toInf.size() > 1) {
			(*it)++; 
		} else {
			continue;
		}
		it++;

		// Need to be a bit clever here. It could be that there are more than 
		// 2 muons that satisfy the pT conditions, that aren't the 2 most 
		// energetic muons, mu1 and mu2 above. We really need to loop through 
		// all pairs to find a suitable pair that meet the criteria
		// That's what the testMuons function does.

		if (DEBUG) cout << "Testing if SS muons" << endl;
		if (testResults(testMuons(muons17toInf, muons10to17, &checkMuonsPTSS), *it)) {
			it++;
		} else {
			continue;
		}

		if (DEBUG) cout << "Testing if eta OK" << endl;
		if (testResults(testMuons(muons17toInf, muons10to17, &checkMuonsPTSSEta), *it)) {
			it++;
		} else {
			continue;
		}

		if (DEBUG) cout << "Testing if dR OK" << endl;
		std::pair<Track*, Track*> p = testMuons(muons17toInf,
					  muons10to17,
					  &checkMuonsSignal);
		if (testResults(p, *it)) {
			it++;
		} else {
			continue;
		}

		Track* mu1 = p.first;
		Track* mu2 = p.second;

		TLorentzVector mu1Mom, mu2Mom;
		mu1Mom = mu1->P4();
		mu2Mom = mu2->P4();

		/////////////////////////////////
		// Look at tracks around muons //
		/////////////////////////////////

		// Vectors of tracks with pT > 1, within dR < 0.5 of respective muons + other cuts
		// so tk1 is the track nearest to muon1, (may or may not be highest pT, depends if random swapping is on)
		std::vector<Track*> tk1_1;
		std::vector<Track*> tk2_1;

		// same but with pT >2.5
		std::vector<Track*> tk1_2p5;
		std::vector<Track*> tk2_2p5;

		// same but with pT > 2.5, OS to muon
		std::vector<Track*> tk1_2p5_OS;
		std::vector<Track*> tk2_2p5_OS;

		Track *candTk(nullptr);
		for(int a = 0; a < branchTracks->GetEntries(); a++){
			candTk = (Track*) branchTracks->At(a);

			if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
				&& (candTk->PT != mu2->PT)
				&& checkTrackLoose(candTk)
			){

				// Store track in suitable vector
				double dR1 = (candTk->P4()).DeltaR(mu1Mom);
				double dR2 = (candTk->P4()).DeltaR(mu2Mom);

				if (dR1 < 0.5){
					tk1_1.push_back(candTk);
				}
				if (dR2 < 0.5){
					tk2_1.push_back(candTk);
				}

				// For 1-prong candidates, need tighter pT and IP cuts
				if (checkTrackPTTight(candTk)
					&& checkTrackIPTight(candTk)){
					// if (dR1 < 0.5){
					// 	tk1_2p5.push_back(candTk);
					// }
					// if (dR2 < 0.5){
					// 	tk2_2p5.push_back(candTk);
					// }
					if ((candTk->Charge) * (mu1->Charge) < 0) {
						if (dR1 < 0.5){
							tk1_2p5_OS.push_back(candTk);
						}
						if (dR2 < 0.5){
							tk2_2p5_OS.push_back(candTk);
						}
					}					
				}

			} // End of track selection
		} // End of track loop

		// Now pT order the track collections for each muon
		// Not needed unless you want the track itself
		// sortTrackVector(tk1_1);
		// sortTrackVector(tk2_1);
		// sortTrackVector(tk1_2p5);
		// sortTrackVector(tk2_2p5);
		// sortTrackVector(tk1_2p5_OS);
		// sortTrackVector(tk2_2p5_OS);

		if (DEBUG) cout << "Testing if only 1 tk with pT > 1 GeV around each muon" << endl;
		if (tk1_1.size() == 1 && tk2_1.size() == 1 ) {
			(*it)++; 
		} else { 
			continue; 
		}
		it++;

		// SIGNAL SELECTION
		if (DEBUG) cout << "Testing if that 1 tk with pT > 2.5 GeV around each muon" << endl;
		if (tk1_1.size() == 1 && tk2_1.size() == 1 
		&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {
			(*it)++; 
		}
		it++;
	} // end of event loop

	// Print out cuts
	for (auto a : cutCount){
		std::cout << a << std::endl;
	}
	
}
