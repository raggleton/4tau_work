#include <vector>
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
	bool swapMuRandomly = true; // if true, fills plots for mu 1 and 2 randomly. Otherwise, does 1 = leading, 2 = subleading
	
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
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_2.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_20.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_3.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_4.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_5.root");
			chain.Add("QCDb_mu_cleanTk/QCDb_mu_6.root");
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
	TH1D *histNTracks1OS       = new TH1D("hNTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 50,0,5);
	TH1D *histNTracks1         = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 50,0,5);
	TH1D *histNTracks2OS       = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 50,0,5);
	TH1D *histNTracks2         = new TH1D("hNTracks2" ,"Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 50,0,5);

	TH1D *histNTracksCum1OS    = new TH1D("hNTracksCum1OS" ,"Cumu Number of tracks about mu1, OS,p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
	TH1D *histNTracksCum1      = new TH1D("hNTracksCum1" ,"Cumu Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
	TH1D *histNTracksCum2OS    = new TH1D("hNTracksCum2OS" ,"Cumu Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);
	TH1D *histNTracksCum2      = new TH1D("hNTracksCum2" ,"Cumu Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);

	TH1D *histMu1Pt            = new TH1D("hMu1Pt", "#mu_{1} p_{T}, no selection ;#mu_{1} p_{T}; N_{events}", 50,0,50.);
	TH1D *histMu2Pt            = new TH1D("hMu2Pt", "#mu_{2} p_{T}, no selection;#mu_{2} p_{T}; N_{events}", 50,0,50.);

	TH1D *histNuPt            = new TH1D("hNuPt", "#nu p_{T}, no selection ;#nu p_{T}; N_{events}", 50,0,50.);
	
	TH1D *histMu1PtSel         = new TH1D("hMu1PtSel", "#mu_{1} p_{T}, selection ;#mu_{1} p_{T}; N_{events}", 50,0,50.);
	TH1D *histMu2PtSel         = new TH1D("hMu2PtSel", "#mu_{2} p_{T}, selection;#mu_{2} p_{T}; N_{events}", 50,0,50.);

	TH1D *histNMu              = new TH1D("hNMu", "No. muons;N mu; N_{events}", 5,0,5);

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
		histNMu->Fill(branchGenMuons->GetEntries());

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
			histRand->Fill(randNum);
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
			} else {
				// cout << "Got > 1 prong!" << endl;
			}

			// cout << "Tau1a has "  << tau1aDaughters.size() << endl;
			// cout << "Tau1b has "  << tau1bDaughters.size() << endl;
			// cout << "Tau2a has "  << tau2aDaughters.size() << endl;
			// cout << "Tau2b has "  << tau2bDaughters.size() << endl;
		} // end if(doSignal)

		////////////////////
		// Muon selection //
		////////////////////
		
		if ( (mu1PT > 10.)
		&& (mu2PT > 10.)
		&& ((mu1->Charge) == (mu2->Charge))
		&& (fabs(origMu1->Eta) < 2.1)
		&& (fabs(origMu2->Eta) < 2.1)
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

			int nAroundMu1 (0), nAroundMu2(0);
			int nTk1(0), nTk25(0);

			n1++;
			n2++;

			// cout << "Track mult: " << branchTracks->GetEntries() << endl;
			// histNTk->Fill(branchTracks->GetEntries());

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
					if (candTk->PT > 2.5){
						nTk25++;
						double dR1 = (candTk->P4()).DeltaR(mu1Mom);
						double dR2 = (candTk->P4()).DeltaR(mu2Mom);

						histNTracks1->Fill(dR1);
						histNTracks2->Fill(dR2);

						histTroubleEtaVsPhi1->Fill(fabs(candTk->Eta - mu1Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu1Mom)));
						histTroubleEtaVsPhi2->Fill(fabs(candTk->Eta - mu2Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu2Mom)));

						// Investigate large dR between trk and muon
						if (candTk->PT > 2.5 && ((dR1 > 0.7 && dR1 < 1.2) || (dR2 > 0.7 && dR2 < 1.2))){
							cout << "CandTk pT: " << candTk->PT << " PID " << candTk->PID << " dR1: " << dR1 << " dR2: " << dR2 << " phi: " << candTk->Phi << " eta: " << candTk->Eta << endl;
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

						if ((candTk->Charge) * (mu1->Charge) < 0){ // only need one if statement because SS muons
							histNTracks1OS->Fill(dR1);
							histNTracks2OS->Fill(dR2);
						}

						// Count number of tracks with pT > 1 within a cone of 0.5 about each muon
						if (dR1 < 0.5)
							nAroundMu1++;
						if (dR2 < 0.5)
							nAroundMu2++;
			
					} //end of 2.5 cut
				} // End of track selection
			} // End of track loop

			// histNTk1->Fill(nTk1);
			// histNTk25->Fill(nTk25);


			if (nAroundMu1==1 && nAroundMu2==1){
				nMuPass++;
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

	histNTracks1->Scale(1./n1);
	histNTracks1OS->Scale(1./n1);
	histNTracks2->Scale(1./n2);
	histNTracks2OS->Scale(1./n2);

	histNTracksCum1   = (TH1D*)histNTracks1->Clone();
	histNTracksCum1OS = (TH1D*)histNTracks1OS->Clone();
	histNTracksCum2   = (TH1D*)histNTracks2->Clone();
	histNTracksCum2OS = (TH1D*)histNTracks2OS->Clone();

	for (int i = 1; i <= histNTracks1->GetNbinsX(); i++){
		histNTracksCum1->SetBinContent(i,histNTracksCum1->GetBinContent(i-1) + histNTracks1->GetBinContent(i));
		histNTracksCum1OS->SetBinContent(i,histNTracksCum1OS->GetBinContent(i-1) + histNTracks1OS->GetBinContent(i));
		histNTracksCum2->SetBinContent(i,histNTracksCum2->GetBinContent(i-1) + histNTracks2->GetBinContent(i));
		histNTracksCum2OS->SetBinContent(i,histNTracksCum2OS->GetBinContent(i-1) + histNTracks2OS->GetBinContent(i));
	}

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
	
	app += "_samePtEta";

	histNMu->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NMu_clean"+app+".pdf").c_str());

	histMu1Pt->Draw("HISTE");
	c.SaveAs((name+"cleanTk/Mu1Pt_clean"+app+".pdf").c_str());
	histMu2Pt->Draw("HISTE");
	c.SaveAs((name+"cleanTk/Mu2Pt_clean"+app+".pdf").c_str());

	histNuPt->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NuPt_clean"+app+".pdf").c_str());
	
	histMu1PtSel->Draw("HISTE");
	c.SaveAs((name+"cleanTk/Mu1PtSel_clean"+app+".pdf").c_str());
	histMu2PtSel->Draw("HISTE");
	c.SaveAs((name+"cleanTk/Mu2PtSel_clean"+app+".pdf").c_str());

	histNTracks1->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks1_NS_clean"+app+".pdf").c_str());
	histNTracks2->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks2_NS_clean"+app+".pdf").c_str());

	histNTracks1OS->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks1_OS_clean"+app+".pdf").c_str());
	histNTracks2OS->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks2_OS_clean"+app+".pdf").c_str());

	histNTracksCum1->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks1Cum_NS_clean"+app+".pdf").c_str());
	histNTracksCum2->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks2Cum_NS_clean"+app+".pdf").c_str());

	histNTracksCum1OS->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks1Cum_OS_clean"+app+".pdf").c_str());
	histNTracksCum2OS->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTracks2Cum_OS_clean"+app+".pdf").c_str());

	histDRMuMu->Draw("HISTE");
	c.SaveAs((name+"cleanTk/DRMuMu_clean"+app+".pdf").c_str());
	histDEtaVsDPhiMuMu->Draw("COLZ");
	c.SaveAs((name+"cleanTk/DEtaVsDPhiMuMu_clean"+app+".pdf").c_str());

	histNTk->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTk_clean"+app+".pdf").c_str());
	histNTk1->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTk1_clean"+app+".pdf").c_str());
	histNTk25->Draw("HISTE");
	c.SaveAs((name+"cleanTk/NTk25_clean"+app+".pdf").c_str());

	if (doSignal){
		histDRa1->Draw("HISTE");
		c.SaveAs((name+"cleanTk/DRa1_clean"+app+".pdf").c_str());
		histDRa2->Draw("HISTE");
		c.SaveAs((name+"cleanTk/DRa2_clean"+app+".pdf").c_str());
		histPID->Draw("HISTE");
		c.SaveAs((name+"cleanTk/PID_clean"+app+".pdf").c_str());
	}

	histTroublePt->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleTkPt_clean"+app+".pdf").c_str());
	histTroublePID->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleTkPID_clean"+app+".pdf").c_str());
	histTroubleEta->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleTkEta_clean"+app+".pdf").c_str());
	histTroublePhi->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleTkPhi_clean"+app+".pdf").c_str());
	histTroubleMatch->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleMatch_clean"+app+".pdf").c_str());
	histTroubleDRMuMu->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleDRMuMu_clean"+app+".pdf").c_str());

	histTroubleDPhiMuMu->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleDPhiMuMu_clean"+app+".pdf").c_str());
	histTroubleDEtaMuMu->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleDEtaMuMu_clean"+app+".pdf").c_str());
	histTroubleMu1Pt->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleMu1Pt_clean"+app+".pdf").c_str());
	histTroubleMu2Pt->Draw("HISTE");
	c.SaveAs((name+"cleanTk/TroubleMu2Pt_clean"+app+".pdf").c_str());

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
	c.SaveAs((name+"cleanTk/TroubleEtaVsPhi1_clean"+app+".pdf").c_str());
	histTroubleEtaVsPhi2->Draw("COLZ");
	problemRing1.Draw("only");
	problemRing2.Draw("only");
	c.SaveAs((name+"cleanTk/TroubleEtaVsPhi2_clean"+app+".pdf").c_str());

	histRand->Draw("HISTE");
	c.SaveAs((name+"cleanTk/RandTest"+app+".pdf").c_str());

	TFile* outFile = TFile::Open((name+"cleanTk/output"+app+".root").c_str(),"RECREATE");

	histNMu->Write("",TObject::kOverwrite);
	histMu1Pt->Write("",TObject::kOverwrite);
	histMu2Pt->Write("",TObject::kOverwrite);
	histMu1PtSel->Write("",TObject::kOverwrite);
	histMu2PtSel->Write("",TObject::kOverwrite);
	histNTracks1->Write("",TObject::kOverwrite);
	histNTracks2->Write("",TObject::kOverwrite);
	histNTracks1OS->Write("",TObject::kOverwrite);
	histNTracks2OS->Write("",TObject::kOverwrite);
	histNTracksCum1->Write("",TObject::kOverwrite);
	histNTracksCum2->Write("",TObject::kOverwrite);
	histNTracksCum1OS->Write("",TObject::kOverwrite);
	histNTracksCum2OS->Write("",TObject::kOverwrite);
	histDRMuMu->Write("",TObject::kOverwrite);
	histNTk->Write("",TObject::kOverwrite);
	histNTk1->Write("",TObject::kOverwrite);
	histNTk25->Write("",TObject::kOverwrite);
	if (doSignal){
		histDRa1->Write("",TObject::kOverwrite);
		histDRa2->Write("",TObject::kOverwrite);
		histPID->Write("",TObject::kOverwrite);
	}

	outFile->Close();

}
