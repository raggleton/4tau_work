#include "commonFunctions.h"

using std::cout;
using std::endl;

void massPlots()
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	bool doSignal = false;
	bool doMu = true; // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = true; // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT = true; // for signal MC - require HLT conditions (29K/5*500K evt) or not (500K evt)

	// Create chain of root trees
	TChain chain("Delphes");
	addInputFiles(&chain, doSignal, doMu, doHLT);

	if (swapMuRandomly) cout << "Swapping mu 1<->2 randomly" << endl;
	else cout << "mu1 has higher pT than mu2" << endl;

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Get pointers to branches used in this analysis
	// Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
	// and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
	// TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchTracks   = treeReader->UseBranch("Track");
	TClonesArray *branchGenMuons = treeReader->UseBranch("OnlyGenMuons");
	TClonesArray *branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray *branchAll      = treeReader->UseBranch("AllParticle");

	//////////////////////
	// Book histograms  //
	//////////////////////
	// Plots for testing invariant mass correlation
	double massBins[6]           = {0,1,2,3,4,10};
	TH1D *histM1              = new TH1D("hM1", "Inv. Mass of 1st system, full selection; m(#mu_{1}-tk) [GeV]; N_{events}",10,0,10);
	TH1D *histM2              = new TH1D("hM2", "Inv. Mass of 2st system, full selection; m(#mu_{2}-tk) [GeV]; N_{events}",10,0,10);

	// MC truth - use actual mu-tk pairs from tau
	TH1D *histM1_truth_0to1   = new TH1D("hM1_truth_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_1to2   = new TH1D("hM1_truth_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_2to3   = new TH1D("hM1_truth_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_3toInf = new TH1D("hM1_truth_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	
	// actual dist using signal selection
	TH1D *histM1_0to1         = new TH1D("hM1_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_1to2         = new TH1D("hM1_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_2to3         = new TH1D("hM1_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_3toInf       = new TH1D("hM1_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);

	// actual dist using sideband selection
	TH1D *histM1_side_0to1    = new TH1D("hM1_side_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_side_1to2    = new TH1D("hM1_side_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_side_2to3    = new TH1D("hM1_side_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_side_3toInf  = new TH1D("hM1_side_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);

	int nMu(0);
	int n1(0), n2(0), nMuPass(0);

	// Loop over all events
	cout << "Nevts : " << numberOfEntries <<endl;
	bool stop = false;
	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		if (branchGenMuons->GetEntries() < 2) continue; // skip if <2 muons!

		//////////////////////////////////////////////////////////////////////
		// First, get the two highest pT muons in the event, store their pT //
		// and pointers to the GenParticles                                 //
		//////////////////////////////////////////////////////////////////////
		
		GenParticle *cand(0),*mu1(0), *mu2(0);
		Track *candTk(0);

		double mu1PT(0.), mu2PT(0.);
		for (int i = 0; i < branchGenMuons->GetEntries(); i++){
			cand = (GenParticle*) branchGenMuons->At(i);
			if (cand->PT > mu1PT) {
				mu1 = cand;
				mu1PT = cand->PT;
			}
		}

		for(int j = 0; j < branchGenMuons->GetEntries(); j++){
			cand = (GenParticle*) branchGenMuons->At(j);
			if ((cand->PT > mu2PT) && (cand->PT != mu1->PT)) {
				mu2 = cand;
				mu2PT = cand->PT;
			}
		}

		// Now randomly swap mu1 - mu2
		GenParticle *origMu1(0), *origMu2(0);
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

		//////////////////////////////////////////////////////
		// Get the hard interaction particles for signal MC //
		// No selection cuts applied (only >=2 muons)       //
		//////////////////////////////////////////////////////
		
		GenParticle *charged1a(0);
		GenParticle *charged1b(0);
		GenParticle *charged2a(0);
		GenParticle *charged2b(0);

		if (doSignal) {
			GenParticle *a1(0), *a2(0);
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
			GenParticle *tau1a(0), *tau1b(0), *tau2a(0), *tau2b(0);
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
				GenParticle* muTruth1(0);
				GenParticle* trackTruth1(0);
				GenParticle* muTruth2(0);
				GenParticle* trackTruth2(0);

				// Assign charged products to be mu or track
				bool truth1HasMu = assignMuonAndTrack(muTruth1, trackTruth1, *charged1a, *charged1b);				
				bool truth2HasMu = assignMuonAndTrack(muTruth2, trackTruth2, *charged2a, *charged2b);

				// NOTE: muons are NOT pT ordered

				if (!truth1HasMu || !truth2HasMu) {
					// cout << "Problem, no truth mu for 1 and/or 2!" << endl;
				} else { 
					// Do m1 distribution in bins of m2 - for MC truth (is it actually correlated?)
					double m1(0.);
					double m2(0.);
					
					// Assign m1 to higher pT muon
					if (muTruth1->PT > muTruth2->PT) {
						m1 = (muTruth1->P4()+trackTruth1->P4()).M();
						m2 = (muTruth2->P4()+trackTruth2->P4()).M();
					} else {
						m2 = (muTruth1->P4()+trackTruth1->P4()).M();
						m1 = (muTruth2->P4()+trackTruth2->P4()).M();
					}

					// Randomly swap trk-mu pairs 1<->2 if desired
					if(swapMuRandomly){
						double randNum = (double)rand() / RAND_MAX;
						if (randNum > 0.5){
							double tmp = m2;
							m2 = m1;
							m1 = tmp;
						}
					}

					// plot mu-tk system properties (MC truth)
					// cout << m1 << "     " << m2 << endl;
					if(m2 < 1.)
						histM1_truth_0to1->Fill(m1);
					else if (m2 < 2.)
						histM1_truth_1to2->Fill(m1);
					else if (m2 < 3.)
						histM1_truth_2to3->Fill(m1);
					else
						histM1_truth_3toInf->Fill(m1);
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
		
		if ((mu1PT > 17.)
		&& (mu2PT > 10.)
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

			Track *candTk(0);
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
					}
				} // End of track selection
			} // End of track loop

			// SIGNAL SELECTION
			if (tk1_1.size() == 1 && tk2_1.size() == 1 
			&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {
				
				TLorentzVector track1Mom=tk1_2p5_OS[0]->P4();
				TLorentzVector track2Mom=tk2_2p5_OS[0]->P4();

				double m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
				double m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

				histM1->Fill(m1);
				histM2->Fill(m2);

				if(m2 < 1.)
					histM1_0to1->Fill(m1);
				else if (m2 < 2.)
					histM1_1to2->Fill(m1);
				else if (m2 < 3.)
					histM1_2to3->Fill(m1);
				else
					histM1_3toInf->Fill(m1);
			}

			// SIDEBAND REGION
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

			// ANOTHER SIDEBAND REGION
			// where at least one muon has an additional track with 1< pT < 2.5,
			// within dR < 0.5. No sign requirement.
			if (tk1_2p5.size() == 1 && tk2_2p5.size() == 1 
				&& ((tk1_1to2p5.size() == 1 && tk1_1to2p5.size() == 1) || (tk1_1to2p5.size() == 0 && tk2_1to2p5.size() == 1) || (tk1_1to2p5.size() == 1 && tk2_1to2p5.size() == 0))){
				
				double m1(0), m2(0);				
				if (tk1_1to2p5.size() == 1)
					m1 = (mu1Mom+tk1_2p5[0]->P4()+tk1_1to2p5[0]->P4()).M();
				else
					m1 = (mu1Mom+tk1_2p5[0]->P4()).M();
				
				if(tk2_1to2p5.size() == 1)
					m2 = (mu2Mom+tk2_2p5[0]->P4()+tk2_1to2p5[0]->P4()).M();
				else
					m2 = (mu2Mom+tk2_2p5[0]->P4()).M();
				
				if(m2 < 1.)
					histM1_side_0to1->Fill(m1);
				else if (m2 < 2.)
					histM1_side_1to2->Fill(m1);
				else if (m2 < 3.)
					histM1_side_2to3->Fill(m1);
				else
					histM1_side_3toInf->Fill(m1);
			}

		} // end of muon selection
		
	} // end of event loop

	TCanvas c;
	std::string app("");
	if (doSignal) {
		app = "sig";
	} else {
		app = "bg";
	}
	if (swapMuRandomly)
		app += "_muRand";
	
	if (doHLT)
		app += "_HLT";
	else
		app += "_NoHLT";

	// Get directory that input file was in - put plots in there
	std::string fullpath = chain.GetFile()->GetDirectory("")->GetName();
	std::vector<std::string> elems;
	split(fullpath, '/', elems);
	std::string directory = elems[0];

	// Get Delphes file config used - last part of directory name
	std::vector<std::string> elems2;
	split(directory, '_', elems2);
	std::string delph = elems2[elems2.size()-1];

	// app += "_samePtEta";

	// Mass plots
	drawHistAndSave(histM1, "HISTE", "M1", directory, app);
	drawHistAndSave(histM2, "HISTE", "M2", directory, app);

	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - signal region;m(#mu_{1}-tk) [GeV]; A.U.", histM1_0to1, histM1_1to2, histM1_2to3, histM1_3toInf, "M1_M2", directory, app);

	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - sideband region (at least 1 #mu has add. tk with p_{T} = (1,2.5));m(#mu_{1}-tk) [GeV]; A.U.", histM1_side_0to1, histM1_side_1to2, histM1_side_2to3, histM1_side_3toInf, "M1_M2_side", directory, app);

	if(doSignal){
		drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - MC truth;m(#mu_{1}-tk) [GeV]; A.U.", histM1_truth_0to1, histM1_truth_1to2, histM1_truth_2to3, histM1_truth_3toInf, "M1_M2_truth", directory, app);
	}

	// TFile* outFile = TFile::Open((name+delph+"/output"+app+".root").c_str(),"RECREATE");

	// histNMu->Write("",TObject::kOverwrite);
	// histMu1Pt->Write("",TObject::kOverwrite);
	// histMu2Pt->Write("",TObject::kOverwrite);
	// histMu1PtSel->Write("",TObject::kOverwrite);
	// histMu2PtSel->Write("",TObject::kOverwrite);
	// histNTracks1->Write("",TObject::kOverwrite);
	// histNTracks2->Write("",TObject::kOverwrite);
	// histNTracks1OS->Write("",TObject::kOverwrite);
	// histNTracks2OS->Write("",TObject::kOverwrite);
	// // histNTracksCum1->Write("",TObject::kOverwrite);
	// // histNTracksCum2->Write("",TObject::kOverwrite);
	// // histNTracksCum1OS->Write("",TObject::kOverwrite);
	// // histNTracksCum2OS->Write("",TObject::kOverwrite);
	// histDRMuMu->Write("",TObject::kOverwrite);
	// histNTk->Write("",TObject::kOverwrite);
	// histNTk1->Write("",TObject::kOverwrite);
	// histNTk25->Write("",TObject::kOverwrite);
	// if (doSignal){
	// 	histDRa1->Write("",TObject::kOverwrite);
	// 	histDRa2->Write("",TObject::kOverwrite);
	// 	histPID->Write("",TObject::kOverwrite);
	// }

	// outFile->Close();

	delete treeReader;
}
