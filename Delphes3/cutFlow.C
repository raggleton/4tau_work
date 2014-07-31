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

	// Hold cut counts and names
	std::vector< std::pair<int, std::string> > cutCount;
	cutCount.push_back(std::make_pair(0,"2+ muons"));
	cutCount.push_back(std::make_pair(0,"1+ mu with pT > 17 GeV"));
	cutCount.push_back(std::make_pair(0,"at least one other mu with pT > 10 GeV"));
	cutCount.push_back(std::make_pair(0,"SS muons"));
	cutCount.push_back(std::make_pair(0,"Mu be within tracker"));
	cutCount.push_back(std::make_pair(0,"dR(mu-mu) > 2"));
	cutCount.push_back(std::make_pair(0,"1 track with pT > 1 GeV within dR of each mu"));
	cutCount.push_back(std::make_pair(0,"Track pT > 2.5 GeV"));
	std::vector<std::pair<int, std::string> >::iterator it = cutCount.begin();

	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		if (DEBUG) cout << "*** Event" << endl;
		// Important! Go back to start of vector
		it = cutCount.begin();

		if (DEBUG) cout << (*it).second << endl;
		if (branchGenMuons->GetEntries() >= 2) {
			(*it).first++; 
		} else continue; // skip if <2 muons!
		it++;

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
		
		// Check mu with pT > 17
		if (DEBUG) cout << (*it).second << endl;
		if (muons17toInf.size() > 0) {
			(*it).first++;
		} else {
			continue;
		}
		it++;

		// Check another mu wiht pT > 10
		if (DEBUG) cout << (*it).second << endl;
		if (muons10to17.size() > 0 || muons17toInf.size() > 1) {
			(*it).first++; 
		} else {
			continue;
		}
		it++;

		// Need to be a bit clever here. It could be that there are more than 
		// 2 muons that satisfy the pT conditions, that aren't the 2 most 
		// energetic muons, mu1 and mu2 above. We really need to loop through 
		// all pairs to find a suitable pair that meet the criteria
		// That's what the testMuons function does.

		// Muon pT + SS
		if (DEBUG) cout << (*it).second << endl;
		if (testResults(testMuons(muons17toInf, muons10to17, &checkMuonsPTSS), (*it).first)) {
			it++;
		} else {
			continue;
		}

		// Muon eta
		if (DEBUG) cout << (*it).second << endl;
		if (testResults(testMuons(muons17toInf, muons10to17, &checkMuonsPTSSEta), (*it).first)) {
			it++;
		} else {
			continue;
		}

		// Muon dR cut
		if (DEBUG) cout << (*it).second << endl;
		std::pair<Track*, Track*> p = testMuons(muons17toInf, muons10to17, &checkMuonsSignal);
		if (testResults(p, (*it).first)) {
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
				fillTrackVectors(candTk, mu1, mu2, &tk1_1, &tk2_1);

				// For 1-prong candidates, need tighter pT and IP cuts
				if (checkTrackTight(candTk)){
					// fillTrackVectors(candTk, mu1, mu2, &tk1_2p5, &tk2_2p5);
					
					// Check tk-mu OS. Only test mu1, as mu1 & 2 are SS
					if (checkTkMuOS(candTk, mu1)) {
						fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS, &tk2_2p5_OS);
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

		if (DEBUG) cout << (*it).second << endl;
		if (tk1_1.size() == 1 && tk2_1.size() == 1 ) {
			(*it).first++; 
		} else { 
			continue; 
		}
		it++;

		// SIGNAL SELECTION
		if (DEBUG) cout << (*it).second << endl;
		if (tk1_1.size() == 1 && tk2_1.size() == 1 
		&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {
			(*it).first++; 
		}
		it++;
	} // end of event loop

	// Print out cuts
	for (auto a : cutCount){
		std::cout << a.second << ": " << a.first << std::endl;
	}
	
}
