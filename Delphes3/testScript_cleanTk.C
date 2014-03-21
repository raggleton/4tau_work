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
using namespace std;

// For sorting track vectors by pT
// Ideally we'd use the templated methods in classes/SortableObject.h ...
bool sortTracksByPT(Track* a, Track* b){ 
	return (a->PT) > (b->PT); 
}

// From the daughters of 2 taus, decide which is track and which is mu
bool assignMuonAndTrack(GenParticle* &mu, GenParticle* &tk, GenParticle &a, GenParticle &b){
	tk = 0;
	mu = 0;

	bool aIsMu = fabs(a.PID) == 13;
	bool bIsMu = fabs(b.PID) == 13;

	// Now look at both a and b
	if (aIsMu && !bIsMu) {
		mu = &a;
		tk = &b;
		return true;
	}

	if (!aIsMu && bIsMu) {
		mu = &b;
		tk = &a;
		return true;
	}

	// the muon is the mu with the higest pT, the other mu becomes a track
	if (aIsMu && bIsMu){
		if (a.PT > b.PT) {
			mu = &a;
			tk = &b;
			return true;
		} else {
			mu = &b;
			tk = &a;
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

void testScript_cleanTk()
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	bool doSignal = true;
	bool doMu = true; // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = false; // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	
	// Create chain of root trees
	TChain chain("Delphes");
	if (doSignal){
		// chain.Add("GG_H_aa.root");
		// chain.Add("sig_test.root");
		// chain.Add("Signal_cleanTk/signal_clean.root");
		// chain.Add("Signal_1prong_cleanTk/signal_1prong_cleanTk.root");
		chain.Add("Signal_1prong_bare/signal_1prong_bare.root");
		// chain.Add("Signal_3prong_cleanTk/signal_3prong_cleanTk.root");
		cout << "Doing signal" << endl;
	} else {
		if (doMu){
			cout << "Doing QCDb_mu" << endl;
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_1.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_2.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_3.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_4.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_5.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_6.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_7.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_8.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_9.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_10.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_11.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_12.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_13.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_14.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_15.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_16.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_17.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_18.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_19.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_20.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_21.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_22.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_23.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_24.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_25.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_26.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_27.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_28.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_29.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_30.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_31.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_32.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_33.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_34.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_35.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_36.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_37.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_38.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_39.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_40.root");
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
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Get pointers to branches used in this analysis
	// Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
	// and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
	// TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchTracks   = treeReader->UseBranch("Track");
	TClonesArray *branchGenMuons = treeReader->UseBranch("OnlyGenMuons");
	TClonesArray *branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray *branchAll      = treeReader->UseBranch("AllParticle");

	// Book histograms
	TH1D *histNTracks1OS       = new TH1D("hNTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks1         = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks2OS       = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);
	TH1D *histNTracks2         = new TH1D("hNTracks2" ,"Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);

	TH1D *histNTracksCum1OS    = new TH1D("hNTracksCum1OS" ,"Cumu Number of tracks about mu1, OS,p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
	TH1D *histNTracksCum1      = new TH1D("hNTracksCum1" ,"Cumu Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
	TH1D *histNTracksCum2OS    = new TH1D("hNTracksCum2OS" ,"Cumu Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);
	TH1D *histNTracksCum2      = new TH1D("hNTracksCum2" ,"Cumu Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);

	TH1D *histMu1Pt            = new TH1D("hMu1Pt", "#mu_{1} p_{T}, no selection ;#mu_{1} p_{T}; N_{events}", 50,0,50.);
	TH1D *histMu2Pt            = new TH1D("hMu2Pt", "#mu_{2} p_{T}, no selection;#mu_{2} p_{T}; N_{events}", 50,0,50.);

	TH1D *histTrack1Pt         = new TH1D("hTrack1Pt","Track 1 p_{T}, signal selection; track 1 p_{T} [GeV]; N_{events}", 25,0,25.);
	TH1D *histTrack2Pt         = new TH1D("hTrack2Pt","Track 2 p_{T}, signal selection; track 2 p_{T} [GeV]; N_{events}", 25,0,25.);
	
	// for mu+tk systems
	TH1D *histSys1Pt           = new TH1D("hSys1Pt", "System 1 p_{T}, signal selection ;System 1 p_{T}; N_{events}", 50,0,50.);
	TH1D *histSys2Pt           = new TH1D("hSys2Pt", "System 2 p_{T}, signal selection;System 2 p_{T}; N_{events}", 50,0,50.);
	TH1D *histDRSys            = new TH1D("hDRSys", "#Delta R(Sys_{1}-Sys_{2}), signal selection;#Delta R(Sys_{1}-Sys_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiSys    = new TH2D("hDEtaVsDPhiSys","dPhi vs dEta for system 2 wrt system 1 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// for mu+tk systems - MC Truth
	TH1D *histSys1PtTruth           = new TH1D("hSys1PtTruth", "System 1 p_{T}, signal selection, MC truth;System 1 p_{T}; N_{events}", 50,0,50.);
	TH1D *histSys2PtTruth           = new TH1D("hSys2PtTruth", "System 2 p_{T}, signal selection, MC truth;System 2 p_{T}; N_{events}", 50,0,50.);
	TH1D *histDRSysTruth            = new TH1D("hDRSysTruth", "#Delta R(Sys_{1}-Sys_{2}), signal selection, MC truth;#Delta R(Sys_{1}-Sys_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiSysTruth    = new TH2D("hDEtaVsDPhiSysTruth","dPhi vs dEta for system 2 wrt system 1, MC Truth ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	TH1D *histNuPt             = new TH1D("hNuPt", "#nu p_{T}, no selection ;#nu p_{T}; N_{events}", 50,0,50.);
	
	TH1D *histMu1PtSel         = new TH1D("hMu1PtSel", "#mu_{1} p_{T}, muon selection, no tk selection;#mu_{1} p_{T}; N_{events}", 50,0,50.);
	TH1D *histMu2PtSel         = new TH1D("hMu2PtSel", "#mu_{2} p_{T}, muon selection, no tk selection;#mu_{2} p_{T}; N_{events}", 50,0,50.);

	// TH1D *histNMu              = new TH1D("hNMu", "No. muons;N mu; N_{events}", 5,0,5);
	// TH1D *histNMu1             = new TH1D("hNMu1", "No. muons about 1;N mu; N_{events}", 5,0,5);
	// TH1D *histNMu2             = new TH1D("hNMu2", "No. muons about 2;N mu; N_{events}", 5,0,5);

	TH1D *histNTk25            = new TH1D("hNTk25", "No. tracks, p_{T} > 2.5 GeV;N_{tk}; N_{events}", 25,0,50);
	TH1D *histNTk1             = new TH1D("hNTk1", "No. tracks, p_{T} > 1 GeV;N_{tk}; N_{events}", 25,0,100);
	TH1D *histNTk              = new TH1D("hNTk", "No. tracks, p_{T} > 0 GeV;N_{tk}; N_{events}", 25,50,250);

	TH1D *histDRMuMu           = new TH1D("hDRMuMu", "#Delta R(#mu-#mu), muon selection;#Delta R(#mu_{1}-#mu_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiMuMu   = new TH2D("hDEtaVsDPhiMuMu","dPhi vs dEta of selection muons ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// Signal-specific ones
	TH1D *histDRa1             = new TH1D("hDRa1","#Delta R(#tau-#tau) 1st a_{1}, no muon selection;#Delta R(#tau-#tau); N_{events}", 10,0,1.);
	TH1D *histDRa2             = new TH1D("hDRa2","#Delta R(#tau-#tau) 2nd a_{1}, no muon selection;#Delta R(#tau-#tau); N_{events}", 10,0,1.);

	// TH1D *histHPt           = new TH1D("hHPt","h_{1,2} p_{T}, what selections?;h_{1,2} p_{T} [GeV]; N_{events}",25,0,50); // Isn't included in MC file from Calchep :(

	TH1D *histPID              = new TH1D("hPID","PID of tau 1-prong; PID; N_{events}", 350,0,350);

	TH1D *histRand             = new TH1D("hRand","Testing rand(); Value; N_{events}", 100,0,1);
	
	TH1D *histTroublePt        = new TH1D("hTroublePt","p_{T} for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,0,25.);
	TH1D *histTroublePID       = new TH1D("hTroublePID","PDGID for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 2200,0,2200);
	TH1D *histTroubleEta       = new TH1D("hTroubleEta","#eta for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,0,5.);
	TH1D *histTroublePhi       = new TH1D("hTroublePhi","#phi for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,-TMath::Pi(),TMath::Pi());
	TH1D *histTroubleMatch     = new TH1D("hTroubleMatch","Whether trks with #Delta R(#mu1-trk) = 0.7 - 1.2 match to any of tau decay products ; Value; N_{events}", 2,0,2.);
	TH1D *histTroubleDRMuMu    = new TH1D("hTroubleDRMuMu","#Delta #R (#mu-#mu) for trks #Delta R(#mu1-trk) = 0.7 - 1.2  ; Value; N_{events}", 20,0,TMath::Pi());
	TH1D *histTroubleDPhiMuMu  = new TH1D("hTroubleDPhiMuMu","#Delta #phi (#mu-#mu) for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 20,0,TMath::Pi());
	TH1D *histTroubleDEtaMuMu  = new TH1D("hTroubleDEtaMuMu","#Delta #eta (#mu-#mu) for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 20,-5,5);
	TH1D *histTroubleMu1Pt     = new TH1D("hTroubleMu1Pt","Mu1 p_{T} for trks #Delta R(#mu1-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,0,35);
	TH1D *histTroubleMu2Pt     = new TH1D("hTroubleMu2Pt","Mu2 p_{T} for trks #Delta R(#mu2-trk) = 0.7 - 1.2 ; Value; N_{events}", 20,0,25);
	TH2D *histTroubleEtaVsPhi1 = new TH2D("hTroubleEtaVsPhi1","dPhi vs dEta of tracks (>2.5 GeV) vs muon 1 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());
	TH2D *histTroubleEtaVsPhi2 = new TH2D("hTroubleEtaVsPhi2","dPhi vs dEta of tracks (>2.5 GeV) vs muon 2 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// Plots for testing invariant mass correlation
	double massBins[6]         = {0,1,2,3,4,10};
	// MC truth - use actual mu-tk pairs from tau
	TH1D *histM1_truth_0to1    = new TH1D("hM1_truth_0to1","m(tk-#mu_{1}) for m(tk-#mu_{2}) = 0-1 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_1to2    = new TH1D("hM1_truth_1to2","m(tk-#mu_{1}) for m(tk-#mu_{2}) = 1-2 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_2to3    = new TH1D("hM1_truth_2to3","m(tk-#mu_{1}) for m(tk-#mu_{2}) = 2-3 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_3toInf  = new TH1D("hM1_truth_3toInf","m(tk-#mu_{1}) for m(tk-#mu_{2}) > 3 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);
	
	// actual dist using selection
	TH1D *histM1_0to1          = new TH1D("hM1_0to1","m(tk-#mu_{1}) for m(tk-#mu_{2}) = 0-1 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);
	TH1D *histM1_1to2          = new TH1D("hM1_1to2","m(tk-#mu_{1}) for m(tk-#mu_{2}) = 1-2 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);
	TH1D *histM1_2to3          = new TH1D("hM1_2to3","m(tk-#mu_{1}) for m(tk-#mu_{2}) = 2-3 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);
	TH1D *histM1_3toInf        = new TH1D("hM1_3toInf","m(tk-#mu_{1}) for m(tk-#mu_{2}) > 3 GeV; m(tk-#mu_{1}) [GeV]; A.U.",5,massBins);

	int nMu(0);
	int n1(0), n2(0), nMuPass(0);

	// Loop over all events
	cout << "Nevts : " << numberOfEntries <<endl;
	bool stop = false;
	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		// Do at gen particle level.
		// histNMu->Fill(branchGenMuons->GetEntries());

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
			// histRand->Fill(randNum);
			if (randNum > 0.5){
				mu1 = origMu2;
				mu2 = origMu1;
			}
		}

		histMu1Pt->Fill(mu1PT);
		histMu2Pt->Fill(mu2PT);

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
			
				if (fabs(cand->PID)== 12|| fabs(cand->PID)== 14|| fabs(cand->PID)==16 ){
					histNuPt->Fill(cand->PT);
				} 

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

			histDRa1->Fill(tau1aMom.DeltaR(tau1bMom));
			histDRa2->Fill(tau2aMom.DeltaR(tau2bMom));
			
			charged1a = getChargedObject(branchAll, tau1a);
			charged1b = getChargedObject(branchAll, tau1b);
			charged2a = getChargedObject(branchAll, tau2a);
			charged2b = getChargedObject(branchAll, tau2b);
			
			if (charged1a && charged1b && charged2a && charged2b){
				
				histPID->Fill(fabs(charged1a->PID));
				histPID->Fill(fabs(charged1b->PID));
				histPID->Fill(fabs(charged2a->PID));
				histPID->Fill(fabs(charged2b->PID));

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

					// cout << m1 << "     " << m2 << endl;
					if(m2 < 1.)
						histM1_truth_0to1->Fill(m1);
					else if (m2 < 2.)
						histM1_truth_1to2->Fill(m1);
					else if (m2 < 3.)
						histM1_truth_2to3->Fill(m1);
					else
						histM1_truth_3toInf->Fill(m1);

					// plot mu-tk system properties (MC truth)

					// combined mu+tk system
					TLorentzVector sys1 = muTruth1->P4()+trackTruth1->P4();
					TLorentzVector sys2 = muTruth2->P4()+trackTruth2->P4();

					// plot the pT, dR, dEta, DPhi of two systems
					histSys1PtTruth->Fill(sys1.Pt());
					histSys2PtTruth->Fill(sys2.Pt());
					histDRSysTruth->Fill(sys1.DeltaR(sys2));
					histDEtaVsDPhiSysTruth->Fill(fabs(sys1.Eta() - sys2.Eta()),fabs(sys1.DeltaPhi(sys2)));

				}
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
			histDRMuMu->Fill(mu1Mom.DeltaR(mu2Mom));
			histDEtaVsDPhiMuMu->Fill(fabs(origMu1->Eta-origMu2->Eta),fabs(mu1Mom.DeltaPhi(mu2Mom)));

			histMu1PtSel->Fill(mu1PT);
			histMu2PtSel->Fill(mu2PT);

			/////////////////////////////////
			// Look at tracks around muons //
			/////////////////////////////////

			int n1AroundMu1 (0), n1AroundMu2(0); // count tracks with pT > 1
			int n25AroundMu1 (0), n25AroundMu2(0); // count tracks with pT > 2.5
			int nTk1(0), nTk25(0);

			n1++;
			n2++;

			// cout << "Track mult: " << branchTracks->GetEntries() << endl;
			// histNTk->Fill(branchTracks->GetEntries());

			// The two tracks
			Track *track1 = new Track();
			Track *track2 = new Track();

			for(int a = 0; a < branchTracks->GetEntries(); a++){
				candTk = (Track*) branchTracks->At(a);

				if ( (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
					&& (candTk->PT != mu2->PT)
					&& (candTk->PT > 1.)
					&& (fabs(candTk->Z) < 1.) //dz < 1mm
					&& ((pow(candTk->X,2)+pow(candTk->Y,2)) < 1.) //dxy < 1mm
					&& (fabs(candTk->Eta)<3)
				){

					nTk1++;
					double dR1 = (candTk->P4()).DeltaR(mu1Mom);
					double dR2 = (candTk->P4()).DeltaR(mu2Mom);
	
					// Also count the track with pT > 1
					if (dR1 < 0.5){
						n1AroundMu1++;
					}
					if (dR2 < 0.5){
						n1AroundMu2++;
					}

					if (candTk->PT > 2.5){

						nTk25++;

						histNTracks1->Fill(dR1);
						histNTracks2->Fill(dR2);

						histTroubleEtaVsPhi1->Fill(fabs(candTk->Eta - mu1Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu1Mom)));
						histTroubleEtaVsPhi2->Fill(fabs(candTk->Eta - mu2Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu2Mom)));

						if ((candTk->Charge) * (mu1->Charge) < 0){ // only need one if statement because SS muons
							histNTracks1OS->Fill(dR1);
							histNTracks2OS->Fill(dR2);
							if ((dR1 < 0.5) && (candTk->PT > track1->PT))
								track1 = candTk;
							if ((dR2 < 0.5) && (candTk->PT > track2->PT))
								track2 = candTk;
						}

						// Count number of tracks with pT > 2.5 within a cone of 0.5 about each muon
						if (dR1 < 0.5)
							n25AroundMu1++;
						if (dR2 < 0.5)
							n25AroundMu2++;

						// Investigate peaked dR between trk and muon, ~ 1
						// If muons are pT ordered (swapMuRandonly = false) then only want the tk about mu 2
						// If randomly ordered, then want tk about both
						if ((swapMuRandomly && ((dR1 > 0.7 && dR1 < 1.4) || (dR2 > 0.7 && dR2 < 1.4))) || (!swapMuRandomly && (dR2 > 0.7 && dR2 < 1.4)) ){
							// cout << "CandTk pT: " << candTk->PT << " PID " << candTk->PID << " dR1: " << dR1 << " dR2: " << dR2 << " phi: " << candTk->Phi << " eta: " << candTk->Eta << endl;
							histTroublePt->Fill(candTk->PT);
							histTroublePID->Fill(fabs(candTk->PID));
							histTroubleEta->Fill(candTk->Eta);
							histTroublePhi->Fill(candTk->Phi);
							histTroubleDRMuMu->Fill(mu1Mom.DeltaR(mu2Mom));
							histTroubleDPhiMuMu->Fill(mu1Mom.DeltaPhi(mu2Mom));
							histTroubleDEtaMuMu->Fill(fabs(mu1Mom.Eta() - mu2Mom.Eta()));
							histTroubleMu1Pt->Fill(mu1PT);
							histTroubleMu2Pt->Fill(mu2PT);
							if (charged1a && charged1b && charged2a && charged2b){
								if (candTk->PT == charged1a->PT || candTk->PT == charged1b->PT || candTk->PT == charged2a->PT || candTk->PT == charged2b->PT)
									histTroubleMatch->Fill(1);
								else
									histTroubleMatch->Fill(0);
							}
							/*for(int j = 0; j < branchAll->GetEntries(); j++){
								candTrouble = (GenParticle*) branchAll->At(j);
								cout << j << " PID: " << candTrouble->PID << " Mother1: " << candTrouble->M1 << " Mother2: " << candTrouble->M2 << " Daughter 1: " << candTrouble->D1 << " Daughter 2: " << candTrouble->D2;
								if ((candTrouble->PT == candTk->PT) && (candTrouble->Eta == candTk->Eta)) {
									cout << " TROUBLE TRACK <<<<<<<<<<<<<<<<<<<<<<<";
								}
								cout << endl;
							}
							stop = true;*/
						}
					} //end of 2.5 cut
				} // End of track selection
			} // End of track loop

			// histNTk1->Fill(nTk1);
			// histNTk25->Fill(nTk25);

			// signal selection - only 1 track with pT > 1, and that track must have pT > 2.5
			if (n1AroundMu1==1 && n1AroundMu2==1 && n25AroundMu1==1 && n25AroundMu1==1){
				nMuPass++;

				histTrack1Pt->Fill(track1->PT);
				histTrack2Pt->Fill(track2->PT);

				TLorentzVector track1Mom=track1->P4();
				TLorentzVector track2Mom=track2->P4();

				// combined mu+tk system
				TLorentzVector sys1 = mu1Mom+track1Mom;
				TLorentzVector sys2 = mu2Mom+track2Mom;

				// plot the pT, dR, dEta, DPhi of two systems
				histSys1Pt->Fill(sys1.Pt());
				histSys2Pt->Fill(sys2.Pt());
				histDRSys->Fill(sys1.DeltaR(sys2));
				histDEtaVsDPhiSys->Fill(fabs(sys1.Eta() - sys2.Eta()),fabs(sys1.DeltaPhi(sys2)));

				// Do m1 in bins of m2
				double m1 = (mu1Mom+track1Mom).M();
				double m2 = (mu2Mom+track2Mom).M();
				cout << m1 << "     " << m2 << endl;
				if(m2 < 1.)
					histM1_0to1->Fill(m1);
				else if (m2 < 2.)
					histM1_1to2->Fill(m1);
				else if (m2 < 3.)
					histM1_2to3->Fill(m1);
				else
					histM1_3toInf->Fill(m1);
			}

			// sideband region where you have 2+ tracks about mu2
			if (n1AroundMu1==1 && n1AroundMu2==1 && n25AroundMu1==1 && n25AroundMu2 >= 2){

			}

		} // end of muon selection
		
		// clean up memory - should prob use smart pointers here
		// delete cand;
		// delete mu1;
		// delete mu2;
		// delete candTk;
		// delete origMu1;
		// delete origMu2;

	} // end of event loop

	histNTracks1->Scale(1./histNTracks1->Integral());
	histNTracks1OS->Scale(1./histNTracks1OS->Integral());
	histNTracks2->Scale(1./histNTracks2->Integral());
	histNTracks2OS->Scale(1./histNTracks2OS->Integral());

	// histNTracksCum1   = (TH1D*)histNTracks1->Clone();
	// histNTracksCum1OS = (TH1D*)histNTracks1OS->Clone();
	// histNTracksCum2   = (TH1D*)histNTracks2->Clone();
	// histNTracksCum2OS = (TH1D*)histNTracks2OS->Clone();

	// for (int i = 1; i <= histNTracks1->GetNbinsX(); i++){
	// 	histNTracksCum1->SetBinContent(i,histNTracksCum1->GetBinContent(i-1) + histNTracks1->GetBinContent(i));
	// 	histNTracksCum1OS->SetBinContent(i,histNTracksCum1OS->GetBinContent(i-1) + histNTracks1OS->GetBinContent(i));
	// 	histNTracksCum2->SetBinContent(i,histNTracksCum2->GetBinContent(i-1) + histNTracks2->GetBinContent(i));
	// 	histNTracksCum2OS->SetBinContent(i,histNTracksCum2OS->GetBinContent(i-1) + histNTracks2OS->GetBinContent(i));
	// }

	cout << "n1: " << n1 << endl;
	cout << "n2: " << n2 << endl;
	cout << "nMuPass: " << nMuPass << endl;

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
	
	std::string delph="bare"; // which Delphes config was used
	// app += "_samePtEta";

	// histNMu->Draw("HISTE");
	// c.SaveAs((name+delph+"/NMu_"+delph+app+".pdf").c_str());
	// histNMu1->Draw("HISTE");
	// c.SaveAs((name+delph+"/NMu1_"+delph+app+".pdf").c_str());
	// histNMu2->Draw("HISTE");
	// c.SaveAs((name+delph+"/NMu2_"+delph+app+".pdf").c_str());

	histMu1Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/Mu1Pt_"+delph+app+".pdf").c_str());
	histMu2Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/Mu2Pt_"+delph+app+".pdf").c_str());

	histTrack1Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/Track1Pt_"+delph+app+".pdf").c_str());
	histTrack2Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/Track2Pt_"+delph+app+".pdf").c_str());

	histSys1Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/Sys1Pt_"+delph+app+".pdf").c_str());
	histSys2Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/Sys2Pt_"+delph+app+".pdf").c_str());
	histDRSys->Draw("HISTE");
	c.SaveAs((name+delph+"/DRSys_"+delph+app+".pdf").c_str());
	histDEtaVsDPhiSys->Draw("COLZ");
	c.SaveAs((name+delph+"/DEtaVsDPhiSys_"+delph+app+".pdf").c_str());

	histSys1PtTruth->Draw("HISTE");
	c.SaveAs((name+delph+"/Sys1PtTruth_"+delph+app+".pdf").c_str());
	histSys2PtTruth->Draw("HISTE");
	c.SaveAs((name+delph+"/Sys2PtTruth_"+delph+app+".pdf").c_str());
	histDRSysTruth->Draw("HISTE");
	c.SaveAs((name+delph+"/DRSysTruth_"+delph+app+".pdf").c_str());
	histDEtaVsDPhiSysTruth->Draw("COLZ");
	c.SaveAs((name+delph+"/DEtaVsDPhiSysTruth_"+delph+app+".pdf").c_str());

	histNuPt->Draw("HISTE");
	c.SaveAs((name+delph+"/NuPt_"+delph+app+".pdf").c_str());
	
	histMu1PtSel->Draw("HISTE");
	c.SaveAs((name+delph+"/Mu1PtSel_"+delph+app+".pdf").c_str());
	histMu2PtSel->Draw("HISTE");
	c.SaveAs((name+delph+"/Mu2PtSel_"+delph+app+".pdf").c_str());

	histNTracks1->Draw("HISTE");
	c.SaveAs((name+delph+"/NTracks1_NS_"+delph+app+".pdf").c_str());
	histNTracks2->Draw("HISTE");
	c.SaveAs((name+delph+"/NTracks2_NS_"+delph+app+".pdf").c_str());

	histNTracks1OS->Draw("HISTE");
	c.SaveAs((name+delph+"/NTracks1_OS_"+delph+app+".pdf").c_str());
	histNTracks2OS->Draw("HISTE");
	c.SaveAs((name+delph+"/NTracks2_OS_"+delph+app+".pdf").c_str());

	// histNTracksCum1->Draw("HISTE");
	// c.SaveAs((name+delph+"/NTracks1Cum_NS_"+delph+app+".pdf").c_str());
	// histNTracksCum2->Draw("HISTE");
	// c.SaveAs((name+delph+"/NTracks2Cum_NS_"+delph+app+".pdf").c_str());

	// histNTracksCum1OS->Draw("HISTE");
	// c.SaveAs((name+delph+"/NTracks1Cum_OS_"+delph+app+".pdf").c_str());
	// histNTracksCum2OS->Draw("HISTE");
	// c.SaveAs((name+delph+"/NTracks2Cum_OS_"+delph+app+".pdf").c_str());

	histDRMuMu->Draw("HISTE");
	c.SaveAs((name+delph+"/DRMuMu_"+delph+app+".pdf").c_str());
	histDEtaVsDPhiMuMu->Draw("COLZ");
	c.SaveAs((name+delph+"/DEtaVsDPhiMuMu_"+delph+app+".pdf").c_str());

	histNTk->Draw("HISTE");
	c.SaveAs((name+delph+"/NTk_"+delph+app+".pdf").c_str());
	histNTk1->Draw("HISTE");
	c.SaveAs((name+delph+"/NTk1_"+delph+app+".pdf").c_str());
	histNTk25->Draw("HISTE");
	c.SaveAs((name+delph+"/NTk25_"+delph+app+".pdf").c_str());

	if (doSignal){
		histDRa1->Draw("HISTE");
		c.SaveAs((name+delph+"/DRa1_"+delph+app+".pdf").c_str());
		histDRa2->Draw("HISTE");
		c.SaveAs((name+delph+"/DRa2_"+delph+app+".pdf").c_str());
		histPID->Draw("HISTE");
		c.SaveAs((name+delph+"/PID_"+delph+app+".pdf").c_str());
	}

	histTroublePt->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleTkPt_"+delph+app+".pdf").c_str());
	histTroublePID->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleTkPID_"+delph+app+".pdf").c_str());
	histTroubleEta->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleTkEta_"+delph+app+".pdf").c_str());
	histTroublePhi->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleTkPhi_"+delph+app+".pdf").c_str());
	histTroubleMatch->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleMatch_"+delph+app+".pdf").c_str());
	histTroubleDRMuMu->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleDRMuMu_"+delph+app+".pdf").c_str());

	histTroubleDPhiMuMu->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleDPhiMuMu_"+delph+app+".pdf").c_str());
	histTroubleDEtaMuMu->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleDEtaMuMu_"+delph+app+".pdf").c_str());
	histTroubleMu1Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleMu1Pt_"+delph+app+".pdf").c_str());
	histTroubleMu2Pt->Draw("HISTE");
	c.SaveAs((name+delph+"/TroubleMu2Pt_"+delph+app+".pdf").c_str());

	TArc problemRing1(0,0,1,0,90);
	problemRing1.SetLineColor(kRed);
	problemRing1.SetLineWidth(2);
	problemRing1.SetFillStyle(0);

	TArc problemRing2(0,0,2,0,90);
	problemRing2.SetLineColor(kRed);
	problemRing2.SetLineWidth(2);
	problemRing2.SetFillStyle(0);

	histTroubleEtaVsPhi1->Draw("COLZ");
	problemRing1.Draw("only");
	c.SaveAs((name+delph+"/TroubleEtaVsPhi1_"+delph+app+".pdf").c_str());
	histTroubleEtaVsPhi2->Draw("COLZ");
	problemRing1.Draw("only");
	problemRing2.Draw("only");
	c.SaveAs((name+delph+"/TroubleEtaVsPhi2_"+delph+app+".pdf").c_str());

	histRand->Draw("HISTE");
	c.SaveAs((name+delph+"/RandTest"+app+".pdf").c_str());

	histM1_truth_0to1->Scale(1./histM1_truth_0to1->Integral());
	histM1_truth_0to1->SetLineColor(kBlack);
	histM1_truth_0to1->Draw("HISTE");
	
	histM1_truth_1to2->Scale(1./histM1_truth_1to2->Integral());
	histM1_truth_1to2->SetLineColor(kRed);
	histM1_truth_1to2->Draw("HISTESAME");
	
	histM1_truth_2to3->Scale(1./histM1_truth_2to3->Integral());
	histM1_truth_2to3->SetLineColor(kGreen);
	histM1_truth_2to3->Draw("HISTESAME");
	
	histM1_truth_3toInf->Scale(1./histM1_truth_3toInf->Integral());
	histM1_truth_3toInf->SetLineColor(kBlue);
	histM1_truth_3toInf->Draw("HISTESAME");

	TLegend leg(0.7,0.7,0.9,0.9);
	leg.AddEntry(histM1_truth_0to1,"m_{2} = 0-1 GeV","l");
	leg.AddEntry(histM1_truth_1to2,"m_{2} = 1-2 GeV","l");
	leg.AddEntry(histM1_truth_2to3,"m_{2} = 2-3 GeV","l");
	leg.AddEntry(histM1_truth_3toInf,"m_{2} > 3 GeV","l");
	leg.Draw();
	c.SaveAs((name+delph+"/M1_truth_"+delph+app+".pdf").c_str());

	histM1_0to1->Scale(1./histM1_0to1->Integral());
	histM1_0to1->SetLineColor(kBlack);
	histM1_0to1->Draw("HISTE");
	
	histM1_1to2->Scale(1./histM1_1to2->Integral());
	histM1_1to2->SetLineColor(kRed);
	histM1_1to2->Draw("HISTESAME");
	
	histM1_2to3->Scale(1./histM1_2to3->Integral());
	histM1_2to3->SetLineColor(kGreen);
	histM1_2to3->Draw("HISTESAME");
	
	histM1_3toInf->Scale(1./histM1_3toInf->Integral());
	histM1_3toInf->SetLineColor(kBlue);
	histM1_3toInf->Draw("HISTESAME");

	leg.Draw();
	c.SaveAs((name+delph+"/M1_"+delph+app+".pdf").c_str());
	
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

}
