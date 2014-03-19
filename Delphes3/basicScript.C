#include <vector>
#include <algorithm>
// #include <iostream>
// #include <utility>

// #include "TROOT.h"
// #include "TSystem.h"
// #include "TApplication.h"

// #include "TString.h"

// #include "TH2.h"
// #include "TH1.h"
// #include "THStack.h"
// #include "TLegend.h"
// #include "TPaveText.h"
// #include "TClonesArray.h"
// #include "TLorentzVector.h"
// #include "TCanvas.h"
// #include "TArc.h"

// #include "classes/DelphesClasses.h"

// #include "external/ExRootAnalysis/ExRootTreeReader.h"
// #include "external/ExRootAnalysis/ExRootTreeWriter.h"
// #include "external/ExRootAnalysis/ExRootTreeBranch.h"
// #include "external/ExRootAnalysis/ExRootResult.h"
// #include "external/ExRootAnalysis/ExRootUtilities.h"
/*
root -l examples/myScript.C\(\"QCDoutput5.root\"\)

for clean tracks ie efficiency = 1, no smearing
 */
using std::cout;
using std::endl;
// using namespace std;

// For sorting track vectors by pT
// Ideally we'd use the templated methods in classes/SortableObject.h ...
bool sortTracksByPT(Track* a, Track* b){ 
	return (a->PT) > (b->PT); 
}

// From the daughters of 2 taus, decide which is track and which is mu
bool assignMuonAndTrack(GenParticle *mu, GenParticle *tk, GenParticle *a, GenParticle *b){
	tk = 0;
	mu = 0;

	bool aIsMu = fabs(a->PID) == 13;
	bool bIsMu = fabs(b->PID) == 13;

	// Now look at both a and b
	if (aIsMu && !bIsMu) {
		mu = a;
		tk = b;
		return true;
	}

	if (!aIsMu && bIsMu) {
		mu = b;
		tk = a;
		return true;
	}

	// the muon is the mu with the higest pT, the other mu becomes a track
	if (aIsMu && bIsMu){
		if (a->PT > b->PT) {
			mu = a;
			tk = b;
			return true;
		} else {
			mu = b;
			tk = a;
			return true;
		}
	}

	// If it gets to here, then neither a nor b is a muon - trouble!
	return false;
}

std::vector<GenParticle*> getTauDaughters(TClonesArray *branchAll, GenParticle *tau) { // get the 3 correct daughters of the tau

	// cout << "get tau daughters" << endl;

	std::vector<GenParticle*> tauDaughters;

	bool foundProduct = false;
	while(!foundProduct){
		if(tau->D1 == tau->D2) { // tau -> tau
			tau = (GenParticle*) branchAll->At(tau->D1);
		} else if (tau->D2-tau->D1 == 1) { //tau->tau+gamma
			if ( fabs(((GenParticle*) branchAll->At(tau->D1))->PID) == 22) {
				tau = (GenParticle*) branchAll->At(tau->D2);
			} else{
				tau = (GenParticle*) branchAll->At(tau->D1);
			}
		} else {
			for(int i = tau->D1; i <= tau->D2; i++){
				tauDaughters.push_back((GenParticle*) branchAll->At(i));
				// tauDaughters.push_back(i);
			}
			if (tauDaughters.size() ==3){
				// cout << "3 products" << endl;
				foundProduct = true;
				return tauDaughters;
			}
		}
	}
}

GenParticle* getChargedObject(TClonesArray* branchAll, GenParticle* tau) { // from tau decay products, get the final stable products

	std::vector<GenParticle*> history; // hold all unique particles in the decay chain (stores event posiiton number)
	std::vector<GenParticle*> current = getTauDaughters(branchAll, tau);
	std::vector<GenParticle*> next; // holds decay products, not nec. all unique
	GenParticle *prong(0); // holds final charged product
	int nCharged = 0; // count charged particles - shouldn't have more than 1!
	bool foundOne = false;

	while (current.size()>0){ // if current > 0 we haven't exhausted all the particles
		for (unsigned a = 0; a < current.size(); a++){
			// Check 1 - is this already in current?
			// Could probably do more efficiently using the unique function on std::vector
			bool alreadyDone = false;

			for (unsigned b = 0; b < a; b++){
				if ((current[a] == current[b]) && (a != 0)) {
					// cout << "--Found a duplicate" << endl;
					alreadyDone = true;
				}
			}

			// Check 2 - is this already in history?
			if (!alreadyDone){
				// cout << "-not in current" << endl;
				for (unsigned c = 0; c < history.size(); c++){
					if ((current[a] == history[c]) && (c!=0)){
						// cout << "--Found a dupicate in history" << endl;
						alreadyDone = true;
					}
				}
			}

			// Check 3 - is this final state?
			// cout << alreadyDone << endl;
			if (!alreadyDone){
				// cout << "-not in history" << endl;

				// Check if final state. Either status 1, or D1 == D2 == -1
				if (current[a]->Status == 1) {

					// Check if charged
					if (current[a]->Charge != 0) {
						// cout << "FINAL PRONG " << current[a]->PID << endl;
						prong = current[a];
						foundOne = true;
						nCharged++;
					}

					history.push_back(current[a]);
					// cout << "pushed into history" << endl;
				} else {
					// Load its daughters no. into next
					for (int d = current[a]->D1; d <= current[a]->D2; d++) {
						// cout << "pushing" << endl;
						next.push_back((GenParticle*) branchAll->At(d));
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

	// Just in case
	history.clear();
	current.clear();
	next.clear();

	if (!foundOne) 
		cout << "No prongs!!!" << endl;

	if (nCharged == 1){
		// cout << "PRONG: " << prong->PID << endl;
		return prong;
	} else {
		// cout << "Damn, more than 1 prong!!!! Got " << nCharged << endl;
		return NULL;
	}

}

void basicScript()
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	bool doSignal = false;
	bool doMu = true; // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = false; // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	
	// Create chain of root trees
	TChain chain("Delphes");
	if (doSignal){
		// chain.Add("GG_H_aa.root");
		// chain.Add("sig_test.root");
		// chain.Add("Signal_cleanTk/signal_clean.root");
		chain.Add("Signal_1prong_cleanTk/signal_1prong_cleanTk.root");
		// chain.Add("Signal_3prong_cleanTk/signal_3prong_cleanTk.root");
		cout << "Doing signal" << endl;
	} else {
		if (doMu){
			cout << "Doing QCDb_mu" << endl;
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_1.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_2.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_3.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_4.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_5.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_6.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_7.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_8.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_9.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_10.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_11.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_12.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_13.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_14.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_15.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_16.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_17.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_18.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_19.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_20.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_21.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_22.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_23.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_24.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_25.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_26.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_27.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_28.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_29.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_30.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_31.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_32.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_33.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_34.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_35.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_36.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_37.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_38.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_39.root");
			// chain.Add("QCDb_mu_cleanTk/QCDb_mu_40.root");
		} else {
			cout << "Doing QCDb" << endl;
			chain.Add("QCDb_cleanTk/QCDb_10.root");
			chain.Add("QCDb_cleanTk/QCDb_2.root");
			chain.Add("QCDb_cleanTk/QCDb_3.root");
			chain.Add("QCDb_cleanTk/QCDb_4.root");
			chain.Add("QCDb_cleanTk/QCDb_5.root");
			chain.Add("QCDb_cleanTk/QCDb_6.root");
			chain.Add("QCDb_cleanTk/QCDb_7.root");
			chain.Add("QCDb_cleanTk/QCDb_8.root");
			chain.Add("QCDb_cleanTk/QCDb_9.root");
		}
	}

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

	// Get pointers to branches used in this analysis
	// Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
	// and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
	// TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchTracks   = treeReader->UseBranch("Track");
	TClonesArray *branchGenMuons = treeReader->UseBranch("OnlyGenMuons");
	TClonesArray *branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray *branchAll      = treeReader->UseBranch("AllParticle");

	// Book histograms
	// TH1D *histNTracks1OS       = new TH1D("hNTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	// TH1D *histNTracks1         = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	// TH1D *histNTracks2OS       = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);
	// TH1D *histNTracks2         = new TH1D("hNTracks2" ,"Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);

	// Loop over all events
	Long64_t numberOfEntries = treeReader->GetEntries();
	cout << "Nevts : " << numberOfEntries <<endl;
	bool stop = false;

	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		if (branchGenMuons->GetEntries() < 2) continue; // skip if <2 muons!

		//////////////////////////////////////////////////////
		// Get the hard interaction particles for signal MC //
		// No selection cuts applied (only >=2 muons)       //
		//////////////////////////////////////////////////////
		
		if (doSignal) {
			GenParticle *a1(0), *a2(0);
			// Get a0s
			for(int j = 0; j < branchAll->GetEntries(); j++){
				GenParticle *candHa = (GenParticle*) branchAll->At(j);
				// cout << j << " ID: " << candHa->PID << " status: " << candHa->Status << endl;
			
				// if (fabs(candHa->PID)== 12|| fabs(candHa->PID)== 14|| fabs(candHa->PID)==16 ){
				// 	histNuPt->Fill(candHa->PT);
				// } 

				// Which is "1" and "2" is arbitrary here.
				if ((fabs(candHa->PID)==36) && (fabs(candHa->Status)==62)){
					if (a1==0){
						a1=candHa;
					} else {
						a2=candHa;
					}
				}
			}

			// Get the tau daughters from a1 and a2 (no pT ordering)
			GenParticle *tau1a(0), *tau1b(0), *tau2a(0), *tau2b(0);
			tau1a = (GenParticle*) branchAll->At(a1->D1);
			tau1b = (GenParticle*) branchAll->At(a1->D2);
			tau2a = (GenParticle*) branchAll->At(a2->D1);
			tau2b = (GenParticle*) branchAll->At(a2->D2);

			// TLorentzVector tau1aMom,tau1bMom, tau2aMom, tau2bMom;
			// tau1aMom = tau1a->P4();
			// tau1bMom = tau1b->P4();
			// tau2aMom = tau2a->P4();
			// tau2bMom = tau2b->P4();

			// histDRa1->Fill(tau1aMom.DeltaR(tau1bMom));
			// histDRa2->Fill(tau2aMom.DeltaR(tau2bMom));
			
			GenParticle *charged1a = getChargedObject(branchAll, tau1a);
			GenParticle *charged1b = getChargedObject(branchAll, tau1b);
			GenParticle *charged2a = getChargedObject(branchAll, tau2a);
			GenParticle *charged2b = getChargedObject(branchAll, tau2b);
			
			// This selects events where each tau only has 1 charged product...
			// dunno what to do about evts where the tau decays into charged things *including* muon
			if (charged1a && charged1b && charged2a && charged2b){
				
				// To hold mu and tracks from each tau
				GenParticle* muTruth1;
				GenParticle* trackTruth1;
				GenParticle* muTruth2;
				GenParticle* trackTruth2;

				// Assign charged products to be mu or track
				bool truth1HasMu = assignMuonAndTrack(muTruth1, trackTruth1, charged1a, charged1b);				
				bool truth2HasMu = assignMuonAndTrack(muTruth2, trackTruth2, charged2a, charged2b);

				// NOTE: muons are NOT pT ordered

				if (!truth1HasMu || !truth2HasMu) {
					cout << "Problem, no truth mu for 1 and/or 2!" << endl;
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

					cout << m1 << "     " << m2 << endl;
					// if(m2 < 1.)
					// 	histM1_truth_0to1->Fill(m1);
					// else if (m2 < 2.)
					// 	histM1_truth_1to2->Fill(m1);
					// else if (m2 < 3.)
					// 	histM1_truth_2to3->Fill(m1);
					// else
					// 	histM1_truth_3toInf->Fill(m1);
				}
			} 
		} // end if(doSignal)


		//////////////////////////////////////////////////////////////////////
		// Now, do more general MC stuff:                                   //
		// get the two highest pT muons in the event, store their pT        //
		// and pointers to the GenParticles                                 //
		//////////////////////////////////////////////////////////////////////
		
		GenParticle *cand(0),*mu1(0), *mu2(0);

		double mu1PT(0.), mu2PT(0.);
		// Get highest pT muon
		for (int i = 0; i < branchGenMuons->GetEntries(); i++){
			cand = (GenParticle*) branchGenMuons->At(i);
			if (cand->PT > mu1PT) {
				mu1 = cand;
				mu1PT = cand->PT;
			}
		}
		// Get 2nd highest pT muon
		for(int j = 0; j < branchGenMuons->GetEntries(); j++){
			cand = (GenParticle*) branchGenMuons->At(j);
			if ((cand->PT > mu2PT) && (cand->PT != mu1->PT)) {
				mu2 = cand;
				mu2PT = cand->PT;
			}
		}

		// Now randomly swap mu1 - mu2 if desired
		GenParticle *origMu1(0), *origMu2(0);
		origMu1 = mu1;
		origMu2 = mu2;
		if (swapMuRandomly){
			double randNum = (double)rand() / RAND_MAX;
			// histRand->Fill(randNum);
			if (randNum > 0.5){
				mu1 = origMu2;
				mu2 = origMu1;
			}
		}

		TLorentzVector mu1Mom, mu2Mom;
		mu1Mom = mu1->P4();
		mu2Mom = mu2->P4();

		////////////////////
		// Muon selection //
		////////////////////
		
		if ((mu1PT < 17.)
		|| (mu2PT < 10.)
		|| ((mu1->Charge) != (mu2->Charge))
		|| (fabs(origMu1->Eta) > 2.1)
		|| (fabs(origMu2->Eta) > 2.4)
		|| ((mu1Mom.DeltaR(mu2Mom)) < 2.)
		){
			continue;
		}

		/////////////////////////////////
		// Look at tracks around muons //
		/////////////////////////////////

		// Vectors of tracks with pT > 1, within dR < 0.5 of respective muons + other cuts
		// so tk1 is the track nearest to muon1, NOT higher pT track
		std::vector<Track*> tk1;
		std::vector<Track*> tk2;

		Track *candTk(0);
		for(int a = 0; a < branchTracks->GetEntries(); a++){
			candTk = (Track*) branchTracks->At(a);

			if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
				&& (candTk->PT != mu2->PT)
				&& (candTk->PT > 1.)
				&& (fabs(candTk->Z) < 1.) //dz < 1mm
				&& ((pow(candTk->X,2)+pow(candTk->Y,2)) < 1.) //dxy < 1mm
				&& (fabs(candTk->Eta)<3)
			){

				double dR1 = (candTk->P4()).DeltaR(mu1Mom);
				double dR2 = (candTk->P4()).DeltaR(mu2Mom);

				if (dR1 < 0.5){
					tk1.push_back(candTk);
				}
				if (dR2 < 0.5){
					tk2.push_back(candTk);
				}

			} // End of track selection
		} // End of track loop

		// Now pT order the track collections
		std::sort(tk1.begin(), tk1.end(), sortTracksByPT);
		std::sort(tk2.begin(), tk2.end(), sortTracksByPT);

		// Can now deal with signal or sidebang regions

		tk1.clear();
		tk2.clear();
	} // end of event loop


	TCanvas c;
	std::string name("");
	std::string app("");
	if (doSignal) {
		// name = "Signal_";
		name = "Signal_1prong_";
		// name = "Signal_3prong_";
		app = "_sig";
	} else {
		app = "_bg";
		if (doMu)
			name = "QCDb_mu_";
		else
			name = "QCDb_";
	}
	if (swapMuRandomly)
		app += "_muRand";
	
	// app += "_samePtEta";

	// histM1_truth_0to1->Scale(1./histM1_truth_0to1->Integral());
	// histM1_truth_0to1->SetLineColor(kBlack);
	// histM1_truth_0to1->Draw("HISTE");
	
	// histM1_truth_1to2->Scale(1./histM1_truth_1to2->Integral());
	// histM1_truth_1to2->SetLineColor(kRed);
	// histM1_truth_1to2->Draw("HISTESAME");
	
	// histM1_truth_2to3->Scale(1./histM1_truth_2to3->Integral());
	// histM1_truth_2to3->SetLineColor(kGreen);
	// histM1_truth_2to3->Draw("HISTESAME");
	
	// histM1_truth_3toInf->Scale(1./histM1_truth_3toInf->Integral());
	// histM1_truth_3toInf->SetLineColor(kBlue);
	// histM1_truth_3toInf->Draw("HISTESAME");

	// TLegend leg(0.7,0.7,0.9,0.9);
	// leg.AddEntry(histM1_truth_0to1,"m_{2} = 0-1 GeV","l");
	// leg.AddEntry(histM1_truth_1to2,"m_{2} = 1-2 GeV","l");
	// leg.AddEntry(histM1_truth_2to3,"m_{2} = 2-3 GeV","l");
	// leg.AddEntry(histM1_truth_3toInf,"m_{2} > 3 GeV","l");
	// leg.Draw();
	// c.SaveAs((name+"cleanTk/M1_truth_clean"+app+".pdf").c_str());

	// histM1_0to1->Scale(1./histM1_0to1->Integral());
	// histM1_0to1->SetLineColor(kBlack);
	// histM1_0to1->Draw("HISTE");
	
	// histM1_1to2->Scale(1./histM1_1to2->Integral());
	// histM1_1to2->SetLineColor(kRed);
	// histM1_1to2->Draw("HISTESAME");
	
	// histM1_2to3->Scale(1./histM1_2to3->Integral());
	// histM1_2to3->SetLineColor(kGreen);
	// histM1_2to3->Draw("HISTESAME");
	
	// histM1_3toInf->Scale(1./histM1_3toInf->Integral());
	// histM1_3toInf->SetLineColor(kBlue);
	// histM1_3toInf->Draw("HISTESAME");

	// leg.Draw();
	// c.SaveAs((name+"cleanTk/M1_clean"+app+".pdf").c_str());
	
	// TFile* outFile = TFile::Open((name+"cleanTk/output"+app+".root").c_str(),"RECREATE");

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

}
