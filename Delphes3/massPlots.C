#include "commonFunctions.h"

using std::cout;
using std::endl;

// int getMassBin(double m){
// 	std::vector<double> massBins {0,1,2,3,10};
// 	for(unsigned i = 1; i < massBins.size(); i++){
// 		if (m < massBins[i])
// 			return i;
// 	}
// 	return 5;
// }

void massPlots(int argc, char* argv[])
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	ProgramOpts pOpts(argc, argv);

	MCsource source     = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	bool doSignal       = pOpts.getSignal(); // for signal or not
	bool doMu           = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = pOpts.getMuOrdering(); // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT          = pOpts.getHLT(); // whether to use MC that has HLT cuts already applied or not.

	// Create chain of root trees
	TChain chain("Delphes");
	addInputFiles(&chain, &pOpts);
	pOpts.printProgramOptions();

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Get pointers to branches used in this analysis
	// Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
	// and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
	// TClonesArray* branchMuon = treeReader->UseBranch("Muon");
	TClonesArray* branchTracks   = treeReader->UseBranch("Track");
	TClonesArray* branchGenMuons = treeReader->UseBranch("OnlyGenMuons");
	// TClonesArray* branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray* branchAll      = treeReader->UseBranch("AllParticle");

	//////////////////////
	// Book histograms  //
	//////////////////////
	// Plots for testing invariant mass correlation
	std::vector<double> massBins {0,1,2,3,10};

	// ------------------------
	// m1 & m2 1D distributions
	// ------------------------

	TH1D* histM1                       = new TH1D("hM1", "Inv. Mass of 1st system, full selection; m(#mu_{1}-tk) [GeV]; N_{events}",10,0,10);
	TH1D* histM2                       = new TH1D("hM2", "Inv. Mass of 2st system, full selection; m(#mu_{2}-tk) [GeV]; N_{events}",10,0,10);

	TH1D* histM1_side_1to2p5           = new TH1D("hM1_side_1to2p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM2_side_1to2p5           = new TH1D("hM2_side_1to2p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",massBins.size()-1,&massBins[0]);
	
	TH1D* histM1_side_1to1p5           = new TH1D("hM1_side_1to1p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM2_side_1to1p5           = new TH1D("hM2_side_1to1p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",massBins.size()-1,&massBins[0]);

	//------------------
	// m1 in bins of m2
	//------------------
	
	// MC truth - use actual mu-tk pairs from tau
	TH1D* histM1_truth_0to1            = new TH1D("hM1_truth_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_truth_1to2            = new TH1D("hM1_truth_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_truth_2to3            = new TH1D("hM1_truth_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_truth_3toInf          = new TH1D("hM1_truth_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	
	// actual dist using signal selection
	TH1D* histM1_0to1                  = new TH1D("hM1_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_1to2                  = new TH1D("hM1_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_2to3                  = new TH1D("hM1_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_3toInf                = new TH1D("hM1_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);

	// actual dist using sideband selection (soft track pT 1-2.5)
	TH1D* histM1_side_1to2p5_0to1      = new TH1D("hM1_side_1to2p5_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_side_1to2p5_1to2      = new TH1D("hM1_side_1to2p5_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_side_1to2p5_2to3      = new TH1D("hM1_side_1to2p5_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_side_1to2p5_3toInf    = new TH1D("hM1_side_1to2p5_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);

	// actual dist using sideband selection (soft track pT 1-1.5)
	TH1D* histM1_side_1to1p5_0to1      = new TH1D("hM1_side_1to1p5_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_side_1to1p5_1to2      = new TH1D("hM1_side_1to1p5_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_side_1to1p5_2to3      = new TH1D("hM1_side_1to1p5_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);
	TH1D* histM1_side_1to1p5_3toInf    = new TH1D("hM1_side_1to1p5_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",massBins.size()-1,&massBins[0]);

	// --------------------------
	// mu1 pT plots in bins of M2
	// --------------------------

	// mu1 pT plots in bins of m2
	TH1D* histMu1Pt_0to1               = new TH1D("hMu1Pt_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_1to2               = new TH1D("hMu1Pt_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_2to3               = new TH1D("hMu1Pt_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_3toInf             = new TH1D("hMu1Pt_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);

	// mu1 pT plots in bins of m2 - MC truth
	TH1D* histMu1Pt_truth_0to1         = new TH1D("hMu1Pt_truth_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_1to2         = new TH1D("hMu1Pt_truth_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_2to3         = new TH1D("hMu1Pt_truth_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_3toInf       = new TH1D("hMu1Pt_truth_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);

	// mu1 pT plots in bins of m2 - sideband (soft track pT 1-2.5)
	TH1D* histMu1Pt_side_1to2p5_0to1   = new TH1D("hMu1Pt_side_1to2p5_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to2p5_1to2   = new TH1D("hMu1Pt_side_1to2p5_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to2p5_2to3   = new TH1D("hMu1Pt_side_1to2p5_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to2p5_3toInf = new TH1D("hMu1Pt_side_1to2p5_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);

	// mu1 pT plots in bins of m2 - sideband (soft track pT 1-1.5)
	TH1D* histMu1Pt_side_1to1p5_0to1   = new TH1D("hMu1Pt_side_1to1p5_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to1p5_1to2   = new TH1D("hMu1Pt_side_1to1p5_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to1p5_2to3   = new TH1D("hMu1Pt_side_1to1p5_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to1p5_3toInf = new TH1D("hMu1Pt_side_1to1p5_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);

	// -------------------
	// 2D plots of m1 V m2
	// -------------------

	// 2D plots of m1 Vs m2 - signal region
	TH2D* histM1vsM2                   = new TH2D("hM1vsM2","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (signal region);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",massBins.size()-1,&massBins[0],massBins.size()-1,&massBins[0]);
	TH2D* histM1timesM1                = new TH2D("hM1timesM2","m(sideband) #times m(sideband) (signal region);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",massBins.size()-1,&massBins[0],massBins.size()-1,&massBins[0]);

	// 2D plots of m1 Vs m2 - sideband (soft track pT 1-2.5)
	TH2D* histM1vsM2_side_1to2p5       = new TH2D("hM1vsM2_side_1to2p5","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",massBins.size()-1,&massBins[0],massBins.size()-1,&massBins[0]);
	TH2D* histM1timesM1_side_1to2p5    = new TH2D("hM1timesM2_side_1to2p5","m(sideband) #times m(sideband) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",massBins.size()-1,&massBins[0],massBins.size()-1,&massBins[0]);

	// 2D plots of m1 Vs m2 - sideband (soft track pT 1-1.5)
	TH2D* histM1vsM2_side_1to1p5       = new TH2D("hM1vsM2_side_1to1p5","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (soft tk p_{T} = 1-1.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",massBins.size()-1,&massBins[0],massBins.size()-1,&massBins[0]);
	TH2D* histM1timesM1_side_1to1p5    = new TH2D("hM1timesM2_side_1to1p5","m(sideband) #times m(sideband) (soft tk p_{T} = 1-1.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",massBins.size()-1,&massBins[0],massBins.size()-1,&massBins[0]);

	// int nMu(0);
	// int n1(0), n2(0), nMuPass(0);

	///////////////////////////
	// Loop over all events  //
	///////////////////////////
	// numberOfEntries = 50000; // for testing only!
	//-------------------------
	cout << "Nevts : " << numberOfEntries << endl;
	bool stop = false; // used to stop the loop, for debugging/testing
	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		if (branchGenMuons->GetEntries() < 2) continue; // skip if <2 muons!

		//////////////////////////////////////////////////////////////////////
		// Now, get the two highest pT muons in the event, store their pT //
		// and pointers to the GenParticles                                 //
		//////////////////////////////////////////////////////////////////////
		
		GenParticle *cand(nullptr),*mu1(nullptr), *mu2(nullptr);
		// Track *candTk(nullptr);

		double muLeadingPT = 0.;
		double muSubLeadingPT = 0.; 
		for (int i = 0; i < branchGenMuons->GetEntries(); i++){
			cand = (GenParticle*) branchGenMuons->At(i);
			if (cand->PT > muLeadingPT) {
				mu1 = cand;
				muLeadingPT = cand->PT;
			}
		}

		for(int j = 0; j < branchGenMuons->GetEntries(); j++){
			cand = (GenParticle*) branchGenMuons->At(j);
			if ((cand->PT > muSubLeadingPT) && (cand->PT != mu1->PT)) {
				mu2 = cand;
				muSubLeadingPT = cand->PT;
			}
		}

		// Now randomly swap mu1 - mu2
		GenParticle *origMu1(nullptr), *origMu2(nullptr);
		origMu1 = mu1;
		origMu2 = mu2;
		if (swapMuRandomly){
			double randNum = (double)rand() / RAND_MAX;
			if (randNum > 0.5){
				mu1 = origMu2;
				mu2 = origMu1;
			}
		}


		TLorentzVector mu1Mom, mu2Mom;
		mu1Mom = mu1->P4();
		mu2Mom = mu2->P4();

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
				cand = (GenParticle*) branchAll->At(j);
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
					if(m2 < 1.){
						histM1_truth_0to1->Fill(m1);
						histMu1Pt_truth_0to1->Fill(muTruth1->PT);
					} else if (m2 < 2.){
						histM1_truth_1to2->Fill(m1);
						histMu1Pt_truth_1to2->Fill(muTruth1->PT);
					} else if (m2 < 3.){
						histM1_truth_2to3->Fill(m1);
						histMu1Pt_truth_2to3->Fill(muTruth1->PT);
					} else{
						histM1_truth_3toInf->Fill(m1);
						histMu1Pt_truth_3toInf->Fill(muTruth1->PT);
					}
				}

				if (!muTruth1) delete muTruth1;
				if (!trackTruth1) delete trackTruth1;
				if (!muTruth2) delete muTruth2;
				if (!trackTruth2) delete trackTruth2;
			} // end if(charged1a...) 
		} // end if(doSignal)

		////////////////////
		// Muon selection //
		////////////////////
		
		if ((muLeadingPT > 17.)
		&& (muSubLeadingPT > 10.)
		&& ((mu1->Charge) == (mu2->Charge))
		&& (fabs(origMu1->Eta) < 2.1)
		// && (fabs(origMu2->Eta) < 2.1)
		&& (fabs(origMu2->Eta) < 2.4)
		&& ((mu1Mom.DeltaR(mu2Mom)) > 2.)
		){
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
			
			// same but with 1 < pT < 2.5 (for sideband)
			std::vector<Track*> tk1_1to2p5;
			std::vector<Track*> tk2_1to2p5;

			// same but with 1 < pT < 1.5 (for sideband)
			std::vector<Track*> tk1_1to1p5;
			std::vector<Track*> tk2_1to1p5;

			Track *candTk(nullptr);
			for(int a = 0; a < branchTracks->GetEntries(); a++){
				candTk = (Track*) branchTracks->At(a);

				if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
					&& (candTk->PT != mu2->PT)
					&& (candTk->PT > 1.)
					&& (fabs(candTk->Z) < 1.) //dz < 1mm
					&& ((pow(candTk->X,2)+pow(candTk->Y,2)) < 1.) //dxy < 1mm
					&& (fabs(candTk->Eta)<2.4)
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

					if (candTk->PT > 2.5){
						if (dR1 < 0.5){
							tk1_2p5.push_back(candTk);
						}
						if (dR2 < 0.5){
							tk2_2p5.push_back(candTk);
						}
						if ((candTk->Charge) * (mu1->Charge) < 0) {
							if (dR1 < 0.5){
								tk1_2p5_OS.push_back(candTk);
							}
							if (dR2 < 0.5){
								tk2_2p5_OS.push_back(candTk);
							}
						}					
					} else {
						if (dR1 < 0.5){
							tk1_1to2p5.push_back(candTk);
						}
						if (dR2 < 0.5){
							tk2_1to2p5.push_back(candTk);
						}
						if (candTk->PT < 1.5){
							if (dR1 < 0.5){
								tk1_1to1p5.push_back(candTk);
							}
							if (dR2 < 0.5){
								tk2_1to1p5.push_back(candTk);
							}
						}
					}
				} // End of track selection criteria
			} // End of track loop

			/////////////////////////
			// SIGNAL SELECTION    //
			/////////////////////////
			// Only 1 track within ∆R < 0.5 of muon has pT > 1, 
			// and that track must have pT > 2.5, and be oppsite charge to muon
			if (tk1_1.size() == 1 && tk2_1.size() == 1 
			&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {
				
				TLorentzVector track1Mom=tk1_2p5_OS[0]->P4();
				TLorentzVector track2Mom=tk2_2p5_OS[0]->P4();

				// random since mu1Mom andmu2Mom are randomly assigned (if selected at top)
				double m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
				double m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

				histM1->Fill(m1);
				histM2->Fill(m2);
				
				// Fill symmetrically to increase stats
				histM1vsM2->Fill(m1,m2);
				// if (getMassBin(m1) != getMassBin(m2)){
				histM1vsM2->Fill(m2,m1);
				
				if(m2 < 1.){
					histM1_0to1->Fill(m1);
					histMu1Pt_0to1->Fill(mu1->PT);
				} else if (m2 < 2.){
					histM1_1to2->Fill(m1);
					histMu1Pt_1to2->Fill(mu1->PT);
				} else if (m2 < 3.){
					histM1_2to3->Fill(m1);
					histMu1Pt_2to3->Fill(mu1->PT);
				} else {
					histM1_3toInf->Fill(m1);
					histMu1Pt_3toInf->Fill(mu1->PT);
				}

			}

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

			////////////////////////////////////////////////////
			// SIDEBAND REGION - for soft tracks 1 < pT < 2.5 //
			////////////////////////////////////////////////////
			// 
			// For 2D plot of N(m1,m2):
			// Each muon must have exactly 1 OS tk with pT > 2.5 within ∆R < 0.5.
			// And *at least* one muon has 1 additional soft track with 1 < pT < 2.5,
			// within dR < 0.5. No sign requirement on additional soft track.
			//
			// For 1D plot of N(m1), N(m2):
			// Each muon must have exactly 1 OS tk with pT > 2.5 within ∆R < 0.5.
			// And *either* muon can have 0 or 1 additional soft track with 1 < pT < 2.5,
			// within dR < 0.5. No sign requirement on additional soft track.
			// 
			// Both regions share same hard track requirements, 
			// but differ in soft track requirements 
			// (1D plot includes case where both have 0 soft tracks, 2D doesn't)
			if (   tk1_1.size() >= 1 && tk1_1.size() < 2 
				&& tk2_1.size() >= 1 && tk2_1.size() < 2 
				&& tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1 
				&& tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1 
				&& tk1_1to2p5.size() <= 1 && tk2_1to2p5.size() <= 1 
				){
					double m1(0), m2(0);		
					m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
					m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

					// Fill 2D m1 vs m2 plot only if at least one muon has a soft track
					if(!(tk1_1to2p5.size() == 0 && tk2_1to2p5.size() == 0)){ 
						// Fill symmetrically to increase stats
						histM1vsM2_side_1to2p5->Fill(m1,m2);
						// if (getMassBin(m1) != getMassBin(m2)){
						histM1vsM2_side_1to2p5->Fill(m2,m1);
						// }
					}

					// Fill 1D plot for the various m2 bins
					histM1_side_1to2p5->Fill(m1);
					histM2_side_1to2p5->Fill(m2);
					if(m2 < 1.) {
						histM1_side_1to2p5_0to1->Fill(m1);
						histMu1Pt_side_1to2p5_0to1->Fill(mu1->PT);
					} else if (m2 < 2.) {
						histM1_side_1to2p5_1to2->Fill(m1);
						histMu1Pt_side_1to2p5_1to2->Fill(mu1->PT);
					} else if (m2 < 3.) {
						histM1_side_1to2p5_2to3->Fill(m1);
						histMu1Pt_side_1to2p5_2to3->Fill(mu1->PT);
					} else {
						histM1_side_1to2p5_3toInf->Fill(m1);
						histMu1Pt_side_1to2p5_3toInf->Fill(mu1->PT);
					}

			} // end of sideband 1to2p5 

			////////////////////////////////////////////////////
			// SIDEBAND REGION - for soft tracks 1 < pT < 1.5 //
			////////////////////////////////////////////////////
			// 
			// For 2D plot of N(m1,m2):
			// Each muon must have exactly 1 OS tk with pT > 2.5 within ∆R < 0.5.
			// And *at least* one muon has 1 additional soft track with 1 < pT < 1.5,
			// within dR < 0.5. No sign requirement on additional soft track.
			//
			// For 1D plot of N(m1), N(m2):
			// Each muon must have exactly 1 OS tk with pT > 2.5 within ∆R < 0.5.
			// And *either* muon can have 0 or 1 additional soft track with 1 < pT < 1.5,
			// within dR < 0.5. No sign requirement on additional soft track.
			// 
			// Both regions share same hard track requirements, 
			// but differ in soft track requirements 
			// (1D plot includes case where both have 0 soft tracks, 2D doesn't)
			if (   tk1_1.size() >= 1 && tk1_1.size() < 2 
				&& tk2_1.size() >= 1 && tk2_1.size() < 2 
				&& tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1 
				&& tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1 
				&& tk1_1to1p5.size() <= 1 && tk2_1to1p5.size() <= 1 
				){
					double m1(0), m2(0);		
					m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
					m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

					// Fill 2D m1 vs m2 plot only if at least one muon has a soft track
					if(!(tk1_1to1p5.size() == 0 && tk2_1to1p5.size() == 0)){ 
						// Fill symmetrically to increase stats
						histM1vsM2_side_1to1p5->Fill(m1,m2);
						// if (getMassBin(m1) != getMassBin(m2)){
						histM1vsM2_side_1to1p5->Fill(m2,m1);
						// }
					}

					// Fill 1D plot for the various m2 bins
					histM1_side_1to1p5->Fill(m1);
					histM2_side_1to1p5->Fill(m2);
					if(m2 < 1.) {
						histM1_side_1to1p5_0to1->Fill(m1);
						histMu1Pt_side_1to1p5_0to1->Fill(mu1->PT);
					} else if (m2 < 2.) {
						histM1_side_1to1p5_1to2->Fill(m1);
						histMu1Pt_side_1to1p5_1to2->Fill(mu1->PT);
					} else if (m2 < 3.) {
						histM1_side_1to1p5_2to3->Fill(m1);
						histMu1Pt_side_1to1p5_2to3->Fill(mu1->PT);
					} else {
						histM1_side_1to1p5_3toInf->Fill(m1);
						histMu1Pt_side_1to1p5_3toInf->Fill(mu1->PT);
					}

			} // end of sideband 1to2p5 

		} // end of muon selection
		
	} // end of event loop

	// Create sum of m1 and m2 1D hists
	TH1D* histM_side_1to2p5 = new TH1D("hM_side_1to2p5","m(#mu-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu-tk) [GeV];A.U.",massBins.size()-1,&massBins[0]);
	histM_side_1to2p5->Add(histM1_side_1to2p5);
	histM_side_1to2p5->Add(histM2_side_1to2p5);
	
	TH1D* histM_side_1to1p5 = new TH1D("hM_side_1to1p5","m(#mu-tk) in sideband (soft tk p_{T} = 1 - 1.5 GeV);m(#mu-tk) [GeV];A.U.",massBins.size()-1,&massBins[0]);
	histM_side_1to1p5->Add(histM1_side_1to1p5);
	histM_side_1to1p5->Add(histM2_side_1to1p5);

	TH1D* histM = new TH1D("hM", "Inv. Mass of system, full selection; m(#mu-tk) [GeV];A.U.",massBins.size()-1,&massBins[0]);
	histM->Add(histM1);
	histM->Add(histM2);

	// Do some normalizing
	normaliseHist(histM1_side_1to2p5);
	normaliseHist(histM2_side_1to2p5);
	normaliseHist(histM_side_1to2p5);
	normaliseHist(histM1vsM2_side_1to2p5);

	normaliseHist(histM1_side_1to1p5);
	normaliseHist(histM2_side_1to1p5);
	normaliseHist(histM_side_1to1p5);
	normaliseHist(histM1vsM2_side_1to1p5);

	normaliseHist(histM1);
	normaliseHist(histM2);
	normaliseHist(histM);
	normaliseHist(histM1vsM2);
	
	// Do corrections plot by making m1*m1 first, then dividing m1vsm2 by m1*m1
	// Don't need to normalise m1timesm1, as histM1_side already normalised
	for(unsigned a = 1; a <= massBins.size()-1; a++){
		for (unsigned b = 1; b <=massBins.size()-1; b++){
			histM1timesM1_side_1to2p5->SetBinContent(a,b,histM_side_1to2p5->GetBinContent(a)*histM_side_1to2p5->GetBinContent(b));
			histM1timesM1_side_1to2p5->SetBinError(a,b,sqrt(pow(histM_side_1to2p5->GetBinContent(b)*histM_side_1to2p5->GetBinError(a),2)
															+pow(histM_side_1to2p5->GetBinContent(a)*histM_side_1to2p5->GetBinError(b),2)));

			histM1timesM1_side_1to1p5->SetBinContent(a,b,histM_side_1to1p5->GetBinContent(a)*histM_side_1to1p5->GetBinContent(b));
			histM1timesM1_side_1to1p5->SetBinError(a,b,sqrt(pow(histM_side_1to1p5->GetBinContent(b)*histM_side_1to1p5->GetBinError(a),2)
															+pow(histM_side_1to1p5->GetBinContent(a)*histM_side_1to1p5->GetBinError(b),2)));

			histM1timesM1->SetBinContent(a,b,histM->GetBinContent(a)*histM->GetBinContent(b));
			histM1timesM1->SetBinError(a,b,sqrt(pow(histM->GetBinContent(b)*histM->GetBinError(a),2)
												+pow(histM->GetBinContent(a)*histM->GetBinError(b),2)));
		}
	}
	TH2D* histM1vsM2_corrections_side_1to2p5 = (TH2D*)histM1vsM2_side_1to2p5->Clone("hM1vsM2_corrections_side_1to2p5");
	histM1vsM2_corrections_side_1to2p5->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband), (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_corrections_side_1to2p5->Divide(histM1timesM1_side_1to2p5);

	TH2D* histM1vsM2_corrections_side_1to1p5 = (TH2D*)histM1vsM2_side_1to1p5->Clone("hM1vsM2_corrections_side_1to1p5");
	histM1vsM2_corrections_side_1to1p5->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband), (soft tk p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_corrections_side_1to1p5->Divide(histM1timesM1_side_1to1p5);

	TH2D* histM1vsM2_corrections = (TH2D*)histM1vsM2->Clone("hM1vsM2_corrections");
	histM1vsM2_corrections->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_corrections->Divide(histM1timesM1);

	TCanvas c;
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

	// Get directory that input file was in - put plots in there
	std::string directory = getDirectory(chain.GetFile());

	// Get Delphes file config used - last part of directory name
	std::string delph = getDelph(directory);

	// Mass plots
	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - signal region;m(#mu_{1}-tk) [GeV]; A.U.", histM1_0to1, histM1_1to2, histM1_2to3, histM1_3toInf, "M1_M2", directory, app);
	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - sideband region (soft tk with p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV]; A.U.", histM1_side_1to2p5_0to1, histM1_side_1to2p5_1to2, histM1_side_1to2p5_2to3, histM1_side_1to2p5_3toInf, "M1_M2_side_1to2p5_1to1p5", directory, app);
	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - sideband region (soft tk with p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV]; A.U.", histM1_side_1to1p5_0to1, histM1_side_1to1p5_1to2, histM1_side_1to1p5_2to3, histM1_side_1to1p5_3toInf, "M1_M2_side_1to1p5_1to1p5", directory, app);
	
	if(doSignal){
		drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - MC truth;m(#mu_{1}-tk) [GeV]; A.U.", histM1_truth_0to1, histM1_truth_1to2, histM1_truth_2to3, histM1_truth_3toInf, "M1_M2_truth", directory, app);
		drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - MC truth;#mu_{1} p_{T} [GeV]; A.U.", histMu1Pt_truth_0to1, histMu1Pt_truth_1to2, histMu1Pt_truth_2to3, histMu1Pt_truth_3toInf, "Mu1Pt_M2_truth", directory, app);
	}

	drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - signal region;#mu_{1} p_{T} [GeV]; A.U.", histMu1Pt_0to1, histMu1Pt_1to2, histMu1Pt_2to3, histMu1Pt_3toInf, "Mu1Pt_M2", directory, app);
	drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - sideband (soft tk p_{T} = 1 - 2.5 GeV);#mu_{1} p_{T} [GeV]; A.U.", histMu1Pt_side_1to2p5_0to1, histMu1Pt_side_1to2p5_1to2, histMu1Pt_side_1to2p5_2to3, histMu1Pt_side_1to2p5_3toInf, "Mu1Pt_M2_side_1to2p5_1to1p5", directory, app);
	drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - sideband (soft tk p_{T} = 1 - 1.5 GeV);#mu_{1} p_{T} [GeV]; A.U.", histMu1Pt_side_1to1p5_0to1, histMu1Pt_side_1to1p5_1to2, histMu1Pt_side_1to1p5_2to3, histMu1Pt_side_1to1p5_3toInf, "Mu1Pt_M2_side_1to1p5_1to1p5", directory, app);

	drawHistAndSave(histM1_side_1to2p5, "HISTE", "M1_side_1to2p5", directory, app);
	drawHistAndSave(histM2_side_1to2p5, "HISTE", "M2_side_1to2p5", directory, app);
	drawHistAndSave(histM1vsM2_side_1to2p5, "colz","M1vsM2_side_1to2p5", directory, app);
	drawHistAndSave(histM1timesM1_side_1to2p5, "colz","M1timesM1_side_1to2p5", directory, app);
	drawHistAndSave(histM1vsM2_corrections_side_1to2p5, "colzTEXTE","M1vsM2_corrections_side_1to2p5", directory, app);

	drawHistAndSave(histM1_side_1to1p5, "HISTE", "M1_side_1to1p5", directory, app);
	drawHistAndSave(histM2_side_1to1p5, "HISTE", "M2_side_1to1p5", directory, app);
	drawHistAndSave(histM1vsM2_side_1to1p5, "colz","M1vsM2_side_1to1p5", directory, app);
	drawHistAndSave(histM1timesM1_side_1to1p5, "colz","M1timesM1_side_1to1p5", directory, app);
	drawHistAndSave(histM1vsM2_corrections_side_1to1p5, "colzTEXTE","M1vsM2_corrections_side_1to1p5", directory, app);

	drawHistAndSave(histM1, "HISTE", "M1", directory, app);
	drawHistAndSave(histM2, "HISTE", "M2", directory, app);
	drawHistAndSave(histM1vsM2, "colz","M1vsM2", directory, app);
	drawHistAndSave(histM1timesM1, "colz","M1timesM1", directory, app);
	drawHistAndSave(histM1vsM2_corrections, "colzTEXTE","M1vsM2_corrections", directory, app);

	TFile* outFile = TFile::Open((directory+"/output_"+delph+"_"+app+".root").c_str(),"UPDATE");

	// Mass plots
	histM1->Write("",TObject::kOverwrite);
	histM2->Write("",TObject::kOverwrite);
	if (doSignal){
		histM1_truth_0to1->Write("",TObject::kOverwrite);
		histM1_truth_1to2->Write("",TObject::kOverwrite);
		histM1_truth_2to3->Write("",TObject::kOverwrite);
		histM1_truth_3toInf->Write("",TObject::kOverwrite);
	}
	histM1_0to1->Write("",TObject::kOverwrite);
	histM1_1to2->Write("",TObject::kOverwrite);
	histM1_2to3->Write("",TObject::kOverwrite);
	histM1_3toInf->Write("",TObject::kOverwrite);
	
	histM1_side_1to2p5_0to1->Write("",TObject::kOverwrite);
	histM1_side_1to2p5_1to2->Write("",TObject::kOverwrite);
	histM1_side_1to2p5_2to3->Write("",TObject::kOverwrite);
	histM1_side_1to2p5_3toInf->Write("",TObject::kOverwrite);
	
	histM1_side_1to1p5_0to1->Write("",TObject::kOverwrite);
	histM1_side_1to1p5_1to2->Write("",TObject::kOverwrite);
	histM1_side_1to1p5_2to3->Write("",TObject::kOverwrite);
	histM1_side_1to1p5_3toInf->Write("",TObject::kOverwrite);

	histMu1Pt_0to1->Write("",TObject::kOverwrite);
	histMu1Pt_1to2->Write("",TObject::kOverwrite);
	histMu1Pt_2to3->Write("",TObject::kOverwrite);
	histMu1Pt_3toInf->Write("",TObject::kOverwrite);
	
	if (doSignal){
		histMu1Pt_truth_0to1->Write("",TObject::kOverwrite);
		histMu1Pt_truth_1to2->Write("",TObject::kOverwrite);
		histMu1Pt_truth_2to3->Write("",TObject::kOverwrite);
		histMu1Pt_truth_3toInf->Write("",TObject::kOverwrite);
	}

	histMu1Pt_side_1to2p5_0to1->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to2p5_1to2->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to2p5_2to3->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to2p5_3toInf->Write("",TObject::kOverwrite);
	
	histMu1Pt_side_1to1p5_0to1->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to1p5_1to2->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to1p5_2to3->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to1p5_3toInf->Write("",TObject::kOverwrite);
	
	histM1_side_1to2p5->Write("",TObject::kOverwrite);
	histM1vsM2_side_1to2p5->Write("",TObject::kOverwrite);
	histM1vsM2_corrections_side_1to2p5->Write("",TObject::kOverwrite);

	histM1_side_1to1p5->Write("",TObject::kOverwrite);
	histM1vsM2_side_1to1p5->Write("",TObject::kOverwrite);
	histM1vsM2_corrections_side_1to1p5->Write("",TObject::kOverwrite);
	
	histM1vsM2->Write("",TObject::kOverwrite);
	histM1vsM2_corrections->Write("",TObject::kOverwrite);

	outFile->Close();

	delete treeReader;
}
