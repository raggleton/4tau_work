#include "commonFunctions.h"
#include "classes/SortableObject.h"
#include "cuts.h"
// #include "tdrstyle.C"

using std::cout;
using std::endl;


/**
 * NB ALL DISTANCES IN MM, ALL ENERGIES IN GeV
 * SEE DelphesHepMCReader.cc
 */


int countParticle(TClonesArray* branchAll, int PID) {
	int nParticle = 0;
	// cout << "**** EVENT *****" << endl;
	for (int a = 0; a < branchAll->GetEntries(); a++){
		GenParticle* candP = (GenParticle*) branchAll->At(a);
		// cout << candP->PID << endl;
		if (fabs(candP->PID) == fabs(PID)) { 
			nParticle++;
		}
	}
	// cout << "**** ENDE *****" << endl;
	return nParticle;
}

/**
 * makes clone fo hist, appends "suffix" to hist name, 
 * write to currently open file
 */
template <typename T>
void makeCopySave(T* h, std::string suffix="_unnormalised") {
	T* h_clone = (T*) h->Clone((h->GetName()+suffix).c_str());
	h_clone->Write("", TObject::kOverwrite);
}

void reweightingVariables(int argc, char* argv[])
{
	// setTDRStyle();

	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	// Get program options from user and store
	ProgramOpts pOpts(argc, argv);

	MCsource source     = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	bool doSignal       = pOpts.getSignal(); // for signal or not
	// bool doMu        = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = pOpts.getMuOrdering(); // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT          = pOpts.getHLT(); // whether to use MC that has HLT cuts already applied or not.
	// bool DEBUG       = pOpts.getVerbose(); // output debug statments
	double deltaR       = pOpts.getdR(); // dR(mu-mu) value to use

	bool do1to1p5 = false; // for additional sideband studies. Slower?

	// Create chain of root trees
	TChain chain("Delphes");
	addInputFiles(&chain, &pOpts);
	pOpts.printProgramOptions();

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

	// Get pointers to branches used in this analysis
	// Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
	// and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
	// TClonesArray* branchMuon = treeReader->UseBranch("Muon");
	TClonesArray* branchTracks   = treeReader->UseBranch("Track");
	// TClonesArray* branchGenMuons = treeReader->UseBranch("OnlyGenMuons"); // GenParticle object
	TClonesArray* branchGenMuons = treeReader->UseBranch("GenMuon"); // Track object
	// TClonesArray* branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray* branchAll      = treeReader->UseBranch("AllParticle");

	//////////////////////
	// Book histograms  //
	//////////////////////
	// Plots for testing invariant mass correlation
	std::vector<double> massBins {0,1,2,3,10};
	int nBinsX = massBins.size()-1;

	// ------------------------
	// Declare hists
	// ------------------------

	TH1D* histDRmutk_hard = new TH1D("hDRmutk_hard", "#DeltaR(#mu-tk) (hard)",40,0,2);
	TH1D* histDRmutk_soft = new TH1D("hDRmutk_soft", "#DeltaR(#mu-tk) (soft)",40,0,2);

	TH1D* histDRmutk_hard_noniso = new TH1D("hDRmutk_hard_noniso", "#DeltaR(#mu-tk) (hard)",40,0,2);
	TH1D* histDRmutk_soft_noniso = new TH1D("hDRmutk_soft_noniso", "#DeltaR(#mu-tk) (soft)",40,0,2);

	// int nMu(0);
	// int n1(0), n2(0), nMuPass(0);

	///////////////////////
	// Loop over events  //
	///////////////////////
	Long64_t numberOfEntries = getNumberEvents(treeReader, &pOpts);
	cout << "Running over " << numberOfEntries << " events" << endl;

	bool stop = false; // used to stop the loop, for debugging/testing
	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		if (branchGenMuons->GetEntries() < 2) continue; // skip if <2 muons! (alhtough pointless for HLT samples)
	
		////////////////////////////////////////////////////////////
		// Get the hard interaction particles for signal MC truth //
		// No selection cuts applied (only >=2 muons)             //
		////////////////////////////////////////////////////////////
		
		GenParticle *charged1a(nullptr);
		GenParticle *charged1b(nullptr);
		GenParticle *charged2a(nullptr);
		GenParticle *charged2b(nullptr);

		if (doSignal) {
			GenParticle *a1(nullptr), *a2(nullptr);
			// Get a0s
			for(int j = 0; j < branchAll->GetEntries(); j++){
				GenParticle* cand = (GenParticle*) branchAll->At(j);
				// cout << j << " ID: " << cand->PID << " status: " << cand->Status << endl;
			
				if ((fabs(cand->PID)==36) && (fabs(cand->Status)==62)){
					if (a1==0){
						a1=cand;
						// cout << "found first a1 at " << j << endl;
					} else {
						// cout << "found second a1 at " << j << endl;
						a2=cand;
					}
				}
			}

			// Get the tau daughters from a1 and a2
			GenParticle *tau1a(nullptr), *tau1b(nullptr), *tau2a(nullptr), *tau2b(nullptr);
			tau1a = (GenParticle*) branchAll->At(a1->D1);
			tau1b = (GenParticle*) branchAll->At(a1->D2);
			tau2a = (GenParticle*) branchAll->At(a2->D1);
			tau2b = (GenParticle*) branchAll->At(a2->D2);

			TLorentzVector tau1aMom,tau1bMom, tau2aMom, tau2bMom;
			tau1aMom = tau1a->P4();
			tau1bMom = tau1b->P4();
			tau2aMom = tau2a->P4();
			tau2bMom = tau2b->P4();

			charged1a = getChargedObject(branchAll, tau1a);
			charged1b = getChargedObject(branchAll, tau1b);
			charged2a = getChargedObject(branchAll, tau2a);
			charged2b = getChargedObject(branchAll, tau2b);
			
			if (charged1a && charged1b && charged2a && charged2b){
				
				// To hold mu and tracks from each tau
				GenParticle* muTruth1(nullptr);
				GenParticle* trackTruth1(nullptr);
				GenParticle* muTruth2(nullptr);
				GenParticle* trackTruth2(nullptr);

				// Assign charged products to be mu or track
				bool truth1HasMu = assignMuonAndTrack(muTruth1, trackTruth1, *charged1a, *charged1b);				
				bool truth2HasMu = assignMuonAndTrack(muTruth2, trackTruth2, *charged2a, *charged2b);

				// NOTE: muons are NOT pT ordered

				if (!truth1HasMu || !truth2HasMu) {
					// cout << "Problem, no truth mu for 1 and/or 2!" << endl;
				} else { 
					
					// Assign system "1" to higher pT muon
					// Swap obj if necessary
					if (muTruth1->PT < muTruth2->PT) {
						GenParticle* tempMu = muTruth1;
						GenParticle* tempTk = trackTruth1;
						muTruth1 = muTruth2;
						trackTruth1 = trackTruth2;
						muTruth2 = tempMu;
						trackTruth2 = tempTk;
					}
					
					// Randomly swap trk-mu pairs 1<->2 if desired
					if(swapMuRandomly){
						double randNum = (double)rand() / RAND_MAX;
						if (randNum > 0.5){
							GenParticle* tempMu = muTruth1;
							GenParticle* tempTk = trackTruth1;
							muTruth1 = muTruth2;
							trackTruth1 = trackTruth2;
							muTruth2 = tempMu;
							trackTruth2 = tempTk;
						}
					}

					// Do m1 distribution in bins of m2 - for MC truth (is it actually correlated?)
					double m1 = (muTruth1->P4()+trackTruth1->P4()).M();
					double m2 = (muTruth2->P4()+trackTruth2->P4()).M();


					// plot mu-tk system properties (MC truth)
					// cout << m1 << "     " << m2 << endl;
					// if(m2 < 1.){
					// 	histM1_truth_0to1->Fill(m1);
					// 	histMu1Pt_truth_0to1->Fill(muTruth1->PT);
					// } else if (m2 < 2.){
					// 	histM1_truth_1to2->Fill(m1);
					// 	histMu1Pt_truth_1to2->Fill(muTruth1->PT);
					// } else if (m2 < 3.){
					// 	histM1_truth_2to3->Fill(m1);
					// 	histMu1Pt_truth_2to3->Fill(muTruth1->PT);
					// } else{
					// 	histM1_truth_3toInf->Fill(m1);
					// 	histMu1Pt_truth_3toInf->Fill(muTruth1->PT);
					// }
				}

				if (!muTruth1) delete muTruth1;
				if (!trackTruth1) delete trackTruth1;
				if (!muTruth2) delete muTruth2;
				if (!trackTruth2) delete trackTruth2;
			} // end if(charged1a...) 
		} // end if(doSignal)


		/////////////////////////////////////////////////////////////////////////
		// Now, get the two highest pT muons in the event that pass selection, // 
		// store pointers to the Track particles and 4-momenta                 //
		// (Use tracks for muons as store more info about position)
		/////////////////////////////////////////////////////////////////////////
		
		// Track *candTk(nullptr);

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
		if (!(muons17toInf.size() >= 1 && (muons17toInf.size() + muons10to17.size()) >= 2)) continue;

		// Sort both vectors by descending pT
		std::sort(muons17toInf.begin(), muons17toInf.end(), sortByPT<Track>);
		std::sort(muons10to17.begin(), muons10to17.end(), sortByPT<Track>);

		// // For muon pT and eta, we want fewer cuts on muons - just SS, IP, Eta
		// std::pair<Track*, Track*> p_loose = testMuons(muons17toInf, muons10to17, &checkMuonsIPSSEta);
		// if (p_loose.first && p_loose.second) {
		// 	histMu1Pt_fine_IPSS->Fill(p_loose.first->PT);
		// 	histMu1Eta_fine_IPSS->Fill(p_loose.first->Eta);
		// 	histMu2Pt_fine_IPSS->Fill(p_loose.second->PT);
		// 	histMu2Eta_fine_IPSS->Fill(p_loose.second->Eta);
		// }		 

		// // For muon pT and eta, we use standard selection, no tracks
		// std::pair<Track*, Track*> p_loose_DR = testMuons(muons17toInf, muons10to17, &checkMuonsIPSSDREta, deltaR);
		// if (p_loose_DR.first && p_loose_DR.second) {
		// 	histMu1Pt_fine_IPSSDR->Fill(p_loose_DR.first->PT);
		// 	histMu1Eta_fine_IPSSDR->Fill(p_loose_DR.first->Eta);
		// 	histMu2Pt_fine_IPSSDR->Fill(p_loose_DR.second->PT);
		// 	histMu2Eta_fine_IPSSDR->Fill(p_loose_DR.second->Eta);
		// 	HardMuonPtSoftMuonPt_DimuonsH->Fill(p_loose_DR.first->PT);
		// }		 

		// For rest of plots we want our muons to pass signal selection
		// Make pairs, see if they pass all cuts (SS, eta, deltaR, dZ, d0)
		// If they do, store in mu1 and mu2 (mu1 has higher pT)
		std::pair<Track*, Track*> p = testMuons(muons17toInf, muons10to17, &checkMuons, deltaR);
		Track* mu1 = p.first;
		Track* mu2 = p.second;

		if (!(p.first && p.second)) continue;


		// cout << " >>>>> Muon pair details: " << endl;
		// cout << " >>>>> Mu1 pT " << mu1->PT << " charge: " << mu1->Charge << " eta " << mu1->Eta << endl;
		// cout << " >>>>> Mu2 pT " << mu2->PT << " charge: " << mu2->Charge << " eta " << mu2->Eta << endl;
		// cout << " >>>>> deltaR " << mu1->P4().DeltaR(mu2->P4()) << endl;

		// Now randomly swap mu1 - mu2
		Track *origMu1(nullptr), *origMu2(nullptr);
		origMu1 = mu1;
		origMu2 = mu2;
		bool swapped = false;
		if (swapMuRandomly){
			// TLorentzVector mu1MomTmp = mu1->P4();
			// TLorentzVector mu2MomTmp = mu2->P4();
			double randNum = (double)rand() / RAND_MAX;
			if (randNum > 0.5){
				swapped = true; 
				mu1 = origMu2;
				mu2 = origMu1;
			}
		}

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

		// same but with pT >2.5, looser (iso) IP cuts
		std::vector<Track*> tk1_2p5_looseip;
		std::vector<Track*> tk2_2p5_looseip;

		// same but with pT > 2.5, OS to muon, , looser (iso) IP cuts
		std::vector<Track*> tk1_2p5_OS_looseip;
		std::vector<Track*> tk2_2p5_OS_looseip;
		
		// same but with 1 < pT < 2.5 (for sideband) (loose/iso IP cuts)
		std::vector<Track*> tk1_1to2p5;
		std::vector<Track*> tk2_1to2p5;
		
		// same but with 1 < pT < 2.5 (for sideband) with tighter IP cuts (tau candidate IP)
		std::vector<Track*> tk1_1to2p5_taucand;
		std::vector<Track*> tk2_1to2p5_taucand;

		// same but with 1 < pT < 1.5 (for sideband)
		std::vector<Track*> tk1_1to1p5;
		std::vector<Track*> tk2_1to1p5;

		Track *candTk(nullptr);
		for(int a = 0; a < branchTracks->GetEntries(); a++){
			candTk = (Track*) branchTracks->At(a);

			if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
				&& (candTk->PT != mu2->PT)
				&& checkTrackLoose(candTk)
			){

				// Store track in suitable vector
				fillTrackVectors(candTk, mu1, mu2, &tk1_1, &tk2_1);


				// 1 prong candiates, 
				if (checkTrackIPTight(candTk)) {
					// for N_trk = 2or3
					fillTrackVectors(candTk, mu1, mu2, &tk1_1to2p5_taucand, &tk2_1to2p5_taucand);						

					if (checkTrackPTTight(candTk)) {
						// tau candidate must have pT > 2.5 GeV
						fillTrackVectors(candTk, mu1, mu2, &tk1_2p5, & tk2_2p5);

						if (checkTkMuOS(candTk, mu1)) {
							fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS, &tk2_2p5_OS);
						}					
					}
				} else {
					
					
					if (checkTrackPTTight(candTk)) {
						// tau candidate must have pT > 2.5 GeV
						fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_looseip, & tk2_2p5_looseip);

						if (checkTkMuOS(candTk, mu1)) {
							fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS_looseip, &tk2_2p5_OS_looseip);
						}					
					} else {
						fillTrackVectors(candTk, mu1, mu2, &tk1_1to2p5, &tk2_1to2p5);
					}

					if (do1to1p5 && candTk->PT < 1.5){
						fillTrackVectors(candTk, mu1, mu2, &tk1_1to1p5, &tk2_1to1p5);
					}
				}
			} // End of track selection criteria
		} // End of track loop


		// Do some preselection here that's common to most regions of consideration
		// 1 OS track with pT > 2.5 within ∆R < 0.5 of muon
		// if (tk1_2p5_OS.size() >= 1 && tk2_2p5_OS.size() >= 1) {
			// sort track vectors by descending pT
			std::sort(tk1_2p5_OS.begin(), tk1_2p5_OS.end(), sortByPT<Track>);
			std::sort(tk2_2p5_OS.begin(), tk2_2p5_OS.end(), sortByPT<Track>);

			// if (tk2_2p5_OS.size() > 1) {
				// cout << tk2_2p5_OS[0]->PT << " - " << tk2_2p5_OS[1]->PT << endl;
			// }
			// do swapping?
			// TLorentzVector totMom = tk1_2p5_OS[0]->P4() + tk2_2p5_OS[0]->P4() + mu1Mom + mu2Mom;
			// TLorentzVector sys1Mom = tk1_2p5_OS[0]->P4() + mu1Mom;
			// TLorentzVector sys2Mom = tk2_2p5_OS[0]->P4() + mu2Mom;


		// } else { 
		// 	continue; 
		// }


		/////////////////////////
		// SIGNAL SELECTION    //
		/////////////////////////
		// Only 1 track within ∆R < 0.5 of muon has pT > 1, 
		// and that track must have pT > 2.5, and be oppsite charge to muon
		if (tk1_1.size() == 1 && tk2_1.size() == 1 )
		{


			if (tk1_2p5.size() == 1 && tk2_2p5.size() == 1
			&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) 
			{
				
				TLorentzVector track1Mom=tk1_2p5_OS[0]->P4();
				TLorentzVector track2Mom=tk2_2p5_OS[0]->P4();

				// random since mu1Mom andmu2Mom are randomly assigned (if selected at top)
				// double m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
				// double m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

				if (swapped){
					histDRmutk_hard->Fill(track2Mom.DeltaR(mu2Mom));
					histDRmutk_soft->Fill(track1Mom.DeltaR(mu1Mom));
				} else {
					histDRmutk_hard->Fill(track1Mom.DeltaR(mu1Mom));
					histDRmutk_soft->Fill(track2Mom.DeltaR(mu2Mom));
				}

			}
		}

		////////////////////////////////////////////////////
		// SIGNAL SELECTION BUT TAUs HAVE LOOSE IP CUTS //
		////////////////////////////////////////////////////
/*
		if (tk1_1.size() == 1 && tk2_1.size() == 1)
		{

			if (tk1_2p5_OS_looseip.size() == 1 && tk2_2p5_OS_looseip.size() == 1 
				&& tk1_2p5_looseip.size() == 1 && tk2_2p5_looseip.size() == 1)
			{
			
				TLorentzVector track1Mom=tk1_2p5_OS_looseip[0]->P4();
				TLorentzVector track2Mom=tk2_2p5_OS_looseip[0]->P4();

				// random since mu1Mom andmu2Mom are randomly assigned (if selected at top)
				double m1 = (mu1Mom+tk1_2p5_OS_looseip[0]->P4()).M();
				double m2 = (mu2Mom+tk2_2p5_OS_looseip[0]->P4()).M();

			}
		}

*/
		// OLD SIDEBAND REGION
		// one muon has 1 tk > 2.5 (OS), other has 2 or 3 ( 1 tk => "mu1", 2/3 tk => "mu2")
		// so NOT pT ordered (although you could do that by eliminating the else if ... bit)

		// if (tk1_1.size() == 1 && (tk2_1.size() == 2 || tk2_1.size() == 3)
		// 	&& tk1_2p5_OS.size() == 1 && (tk2_2p5.size() == 2 || tk2_2p5.size() == 3)){
			
		// 	// mu1Mom has 1 tk, mu2Mom has 2/3 tks, stay as "mu1" and "mu2"
		// 	double m1 = (mu1Mom+tk1_2p5[0]->P4()).M();
		// 	double m2 = (mu2Mom+tk2_2p5[0]->P4()).M();
		// 	if(m2 < 1.)
		// 		histM1_side_0to1->Fill(m1);
		// 	else if (m2 < 2.)
		// 		histM1_side_1to2->Fill(m1);
		// 	else if (m2 < 3.)
		// 		histM1_side_2to3->Fill(m1);
		// 	else
		// 		histM1_3toInf->Fill(m1);
		// } else if (tk2_1.size() == 1 && (tk1_1.size() == 2 || tk1_1.size() == 3)
		// 	&& tk2_2p5_OS.size() == 1 && (tk1_2p5.size() == 2 || tk1_2p5.size() == 3)){

		// 	// mu1Mom has 2/3 tks, mu2Mom has 1 tks, so m1 uses mu2Mom & v.v.
		// 	double m1 = (mu2Mom+tk2_2p5[0]->P4()).M();
		// 	double m2 = (mu1Mom+tk1_2p5[0]->P4()).M();
		// 	if(m2 < 1.)
		// 		histM1_side_0to1->Fill(m1);
		// 	else if (m2 < 2.)
		// 		histM1_side_1to2->Fill(m1);
		// 	else if (m2 < 3.)
		// 		histM1_side_2to3->Fill(m1);
		// 	else
		// 		histM1_side_3toInf->Fill(m1);
		// }

		//////////////////////////////////////////////////////////////////////
		// SIDEBAND REGION - where mu2 has 1, 2, 3 additional soft tracks with tight IP cuts//
		// and mu1 satisfies signal selection
		//////////////////////////////////////////////////////////////////////
		if (   tk1_1.size() == 1 
			&& tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1 
			&& tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1 
			){
			
			double m1(0);		
			TLorentzVector track1Mom=tk1_2p5_OS[0]->P4();
			TLorentzVector track2Mom=tk2_2p5_OS[0]->P4();

			if ( //tk2_1.size() == 2 
				tk2_1to2p5_taucand.size() == 1
				&& tk2_1to2p5.size() == 1
				) {
					// m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
					

				if (swapped){
					histDRmutk_hard_noniso->Fill(track2Mom.DeltaR(mu2Mom));
					histDRmutk_soft_noniso->Fill(track1Mom.DeltaR(mu1Mom));
				} else {
					histDRmutk_hard_noniso->Fill(track1Mom.DeltaR(mu1Mom));
					histDRmutk_soft_noniso->Fill(track2Mom.DeltaR(mu2Mom));
				}

			}			

			if (// tk2_1.size() == 3 
				tk2_1to2p5_taucand.size() == 2
				&& tk2_1to2p5.size() == 2
				) {
					// m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
					
					if (swapped){
						histDRmutk_hard_noniso->Fill(track2Mom.DeltaR(mu2Mom));
						histDRmutk_soft_noniso->Fill(track1Mom.DeltaR(mu1Mom));
					} else {
						histDRmutk_hard_noniso->Fill(track1Mom.DeltaR(mu1Mom));
						histDRmutk_soft_noniso->Fill(track2Mom.DeltaR(mu2Mom));
					}


			}			

/*			if ( //tk2_1.size() == 4 
				tk2_1to2p5_taucand.size() == 3
				&& tk2_1to2p5.size() == 3
				) {
					m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
				}	*/		
		}


		////////////////////////////////////////////////////
		// SIDEBAND REGION - for soft tracks 1 < pT < 2.5 //
		////////////////////////////////////////////////////
		// 
		// For 2D plot of N(m1,m2):
		// Each muon must have exactly 1 OS tk with pT > 2.5 within ∆R < 0.5.
		// And *at least* one muon has 1 or 2 additional soft tracks with 1 < pT < 2.5,
		// within dR < 0.5. No sign requirement on additional soft track.
		// 
		// This matches Control Region A in the PAS (B in AN)
		/*if (   tk1_1.size() >= 1 && tk1_1.size() <= 3 
			&& tk2_1.size() >= 1 && tk2_1.size() <= 3 
			&& tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1 
			&& tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1 
			&& tk1_1to2p5.size() <= 2 && tk2_1to2p5.size() <= 2 
			&& !(tk1_1to2p5.size() == 0 && tk2_1to2p5.size() == 0)
			){
				double m1(0), m2(0);		
				m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
				m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

				// Fill 1D plot for the various m2 bins
				// histM1_side_1to2p5->Fill(m1);
				// histM2_side_1to2p5->Fill(m2);

		} // end of sideband 1to2p5 

		*/


		////////////////////////////////////////////////////
		// SIDEBAND REGION - for soft tracks 1 < pT < 2.5 //
		// AND WHERE TAU CANDIDATE HAS LOOSER IP CUTS
		////////////////////////////////////////////////////
		// 
		// For 2D plot of N(m1,m2):
		// Each muon must have exactly 1 OS tk with pT > 2.5 within ∆R < 0.5.
		// And *at least* one muon has 1 or 2 additional soft tracks with 1 < pT < 2.5,
		// within dR < 0.5. No sign requirement on additional soft track.
		// 
		// This matches Control Region A in the PAS (B in AN)
		/*if (   tk1_1.size() >= 1 && tk1_1.size() <= 3 
			&& tk2_1.size() >= 1 && tk2_1.size() <= 3 
			&& tk1_2p5_looseip.size() == 1 && tk1_2p5_OS_looseip.size() == 1 
			&& tk2_2p5_looseip.size() == 1 && tk2_2p5_OS_looseip.size() == 1 
			&& tk1_1to2p5.size() <= 2 && tk2_1to2p5.size() <= 2 
			&& !(tk1_1to2p5.size() == 0 && tk2_1to2p5.size() == 0)
			){
				double m1(0), m2(0);		
				m1 = (mu1Mom+tk1_2p5_OS_looseip[0]->P4()).M();
				m2 = (mu2Mom+tk2_2p5_OS_looseip[0]->P4()).M();


		} // end of sideband 1to2p5 
*/
		// SIDEBAND - 1D plot for mu1 
		// mu1 has 1 OS track with pT > 2.5, within ∆R < 0.5. 
		// Also has 0 or 1 soft track, 1 < pT < 2.5, within ∆R < 0.5.
		// No sign requirement.
		// No track requirements on mu2, all it has to do is pass the mu selection above.
		// if(tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1 && tk1_1to2p5.size() <= 1){
		// 	double m1(0);		
		// 	m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
		// 	histM1_side_1to2p5->Fill(m1);
		// }

		// SIDEBAND - 1D plot for mu2
		// mu2 has 1 OS track with pT > 2.5, within ∆R < 0.5. 
		// Also has 0 or 1 soft track, 1 < pT < 2.5, within ∆R < 0.5.
		// No sign requirement.
		// No track requirements on mu1, all it has to do is pass the mu selection above.
		// if(tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1 && tk2_1to2p5.size() <= 1){
		// 	double m2(0);		
		// 	m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();
		// 	histM2_side_1to2p5->Fill(m2);
		// }

		////////////////////////////////////////////////////
		// SIDEBAND REGION - for soft tracks 1 < pT < 1.5 //
		////////////////////////////////////////////////////
		//

	} // end of event loop


	// Print out integrals before we normalise
	printIntegral(histDRmutk_hard);
	printIntegral(histDRmutk_soft);
	printIntegral(histDRmutk_hard_noniso);
	printIntegral(histDRmutk_soft_noniso);

	/////////////////
	// PLOT THINGS //
	/////////////////
	std::string app(""); // text to append on end of plot filenames

	if (doSignal) {
		app = "sig";
	} else {
		app = "bg";
	}
	if (swapMuRandomly){
		app += "_muRand";
	}
	if (doHLT) {
		app += "_HLT";
	} else {
		app += "_NoHLT";
	}

	app += "_dR";
	app += boost::lexical_cast<std::string>(deltaR);

	if (source == test)
		app += "_TEST";

	// Get directory that input file was in - put plots in there
	std::string directory = getDirectory(chain.GetFile());

	// Get Delphes file config used - last part of directory name
	std::string delph = getDelph(directory);

	TFile* outFile = TFile::Open((directory+"/output_reweight_"+delph+"_"+app+".root").c_str(),"UPDATE");
	cout << "Writing to " << outFile->GetName() << endl;
	
	// Do some normalizing, make copies beforehand for combination plots
	// makeCopySave(histM1_side_1to2p5);
	// normaliseHist(histM1_side_1to2p5);
	
	drawHistAndSave(histDRmutk_hard, "HISTE", "DRmutk_hard", directory, app);
	drawHistAndSave(histDRmutk_soft, "HISTE", "DRmutk_soft", directory, app);
	drawHistAndSave(histDRmutk_hard_noniso, "HISTE", "DRmutk_hard_noniso", directory, app);
	drawHistAndSave(histDRmutk_soft_noniso, "HISTE", "DRmutk_soft_noniso", directory, app);
	
	//////////////////////////
	// Write hists to file //
	//////////////////////////

	histDRmutk_hard->Write("",TObject::kOverwrite);
	histDRmutk_soft->Write("",TObject::kOverwrite);
	histDRmutk_hard_noniso->Write("",TObject::kOverwrite);
	histDRmutk_soft_noniso->Write("",TObject::kOverwrite);

	if (doSignal){
	}

	outFile->Close();

	delete treeReader;
}
