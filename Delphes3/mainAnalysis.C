#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include "commonFunctions.h"

using std::cout;
using std::endl;

/**
 * Main analysis script to make plots etc. 
 */
void mainAnalysis()
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	bool doSignal = false;
	bool doMu = true; // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = true; // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT = false; // whether to use MC that has HLT cuts already applied or not.

	// Create chain of root trees
	TChain chain("Delphes");
	addInputFiles(&chain, doSignal, doMu, doHLT);

	if (swapMuRandomly) cout << "Swapping mu 1<->2 randomly" << endl;
	else cout << "mu1 has higher pT than mu2" << endl;

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

	//////////////////////
	// Book histograms //
	//////////////////////
	TH1D *histNTracks1           = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks2           = new TH1D("hNTracks2" ,"Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);
	TH1D *histNTracks1OS         = new TH1D("hNTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks2OS         = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);

	TH1D *histNTracksCum1        = new TH1D("hNTracksCum1" ,"Cumulative Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNTracksCum2        = new TH1D("hNTracksCum2" ,"Cumulative Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);
	TH1D *histNTracksCum1OS      = new TH1D("hNTracksCum1OS" ,"Cumulative Number of tracks about mu1, OS,p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNTracksCum2OS      = new TH1D("hNTracksCum2OS" ,"Cumulative Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);

	TH1D *histMu1Pt              = new TH1D("hMu1Pt", "#mu_{1} p_{T}, no selection ;#mu_{1} p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histMu2Pt              = new TH1D("hMu2Pt", "#mu_{2} p_{T}, no selection;#mu_{2} p_{T} [GeV]; N_{events}", 50,0,50.);

	TH1D *histTrack1Pt           = new TH1D("hTrack1Pt","Track 1 p_{T}, signal selection; track 1 p_{T} [GeV]; N_{events}", 25,0,25.);
	TH1D *histTrack2Pt           = new TH1D("hTrack2Pt","Track 2 p_{T}, signal selection; track 2 p_{T} [GeV]; N_{events}", 25,0,25.);
	
	// for mu+tk systems
	TH1D *histSys1Pt             = new TH1D("hSys1Pt", "System 1 p_{T}, signal selection ;System 1 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histSys2Pt             = new TH1D("hSys2Pt", "System 2 p_{T}, signal selection;System 2 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histDRSys              = new TH1D("hDRSys", "#Delta R(Sys_{1}-Sys_{2}), signal selection;#Delta R(Sys_{1}-Sys_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiSys      = new TH2D("hDEtaVsDPhiSys","dPhi vs dEta for system 2 wrt system 1 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// for mu+tk systems - MC Truth
	TH1D *histSys1PtTruth        = new TH1D("hSys1PtTruth", "System 1 p_{T}, signal selection, MC truth;System 1 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histSys2PtTruth        = new TH1D("hSys2PtTruth", "System 2 p_{T}, signal selection, MC truth;System 2 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histDRSysTruth         = new TH1D("hDRSysTruth", "#Delta R(Sys_{1}-Sys_{2}), signal selection, MC truth;#Delta R(Sys_{1}-Sys_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiSysTruth = new TH2D("hDEtaVsDPhiSysTruth","dPhi vs dEta for system 2 wrt system 1, MC Truth ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	TH1D *histNuPt               = new TH1D("hNuPt", "#nu p_{T}, no selection ;#nu p_{T} [GeV]; N_{events}", 50,0,50.);
	
	TH1D *histMu1PtSel           = new TH1D("hMu1PtSel", "#mu_{1} p_{T}, muon selection, no tk selection;#mu_{1} p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histMu2PtSel           = new TH1D("hMu2PtSel", "#mu_{2} p_{T}, muon selection, no tk selection;#mu_{2} p_{T} [GeV]; N_{events}", 50,0,50.);

	// TH1D *histNMu             = new TH1D("hNMu", "No. muons;N mu; N_{events}", 5,0,5);
	// TH1D *histNMu1            = new TH1D("hNMu1", "No. muons about 1;N mu; N_{events}", 5,0,5);
	// TH1D *histNMu2            = new TH1D("hNMu2", "No. muons about 2;N mu; N_{events}", 5,0,5);

	TH1D *histNTk25              = new TH1D("hNTk25", "No. tracks, p_{T} > 2.5 GeV;N_{tk}; N_{events}", 25,0,50);
	TH1D *histNTk1               = new TH1D("hNTk1", "No. tracks, p_{T} > 1 GeV;N_{tk}; N_{events}", 25,0,100);
	TH1D *histNTk                = new TH1D("hNTk", "No. tracks, p_{T} > 0 GeV;N_{tk}; N_{events}", 25,50,250);

	TH1D *histDRMuMu             = new TH1D("hDRMuMu", "#Delta R(#mu-#mu), muon selection;#Delta R(#mu_{1}-#mu_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiMuMu     = new TH2D("hDEtaVsDPhiMuMu","dPhi vs dEta of selection muons ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// Signal-specific ones
	TH1D *histDRa1               = new TH1D("hDRa1","#Delta R(#tau-#tau) 1st a_{1}, no muon selection;#Delta R(#tau-#tau); N_{events}", 10,0,1.);
	TH1D *histDRa2               = new TH1D("hDRa2","#Delta R(#tau-#tau) 2nd a_{1}, no muon selection;#Delta R(#tau-#tau); N_{events}", 10,0,1.);

	TH1D *histHPt                = new TH1D("hHPt","h_{1,2} p_{T}, no selection ;h_{1,2} p_{T} [GeV]; N_{events}",25,0,50); // Isn't included in MC file from Calchep :(

	TH1D *histPID                = new TH1D("hPID","PID of tau 1-prong; PID; N_{events}", 350,0,350);

	TH1D *histRand               = new TH1D("hRand","Testing rand(); Value; N_{events}", 100,0,1);
	
	TH1D *histTroublePt          = new TH1D("hTroublePt","p_{T} for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,0,25.);
	TH1D *histTroublePID         = new TH1D("hTroublePID","PDGID for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 2200,0,2200);
	TH1D *histTroubleEta         = new TH1D("hTroubleEta","#eta for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,0,5.);
	TH1D *histTroublePhi         = new TH1D("hTroublePhi","#phi for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,-TMath::Pi(),TMath::Pi());
	TH1D *histTroubleMatch       = new TH1D("hTroubleMatch","Whether trks with #Delta R(#mu1-trk) = 0.7 - 1.2 match to any of tau decay products ; Value; N_{events}", 2,0,2.);
	TH1D *histTroubleDRMuMu      = new TH1D("hTroubleDRMuMu","#Delta #R (#mu-#mu) for trks #Delta R(#mu1-trk) = 0.7 - 1.2  ; Value; N_{events}", 20,0,TMath::Pi());
	TH1D *histTroubleDPhiMuMu    = new TH1D("hTroubleDPhiMuMu","#Delta #phi (#mu-#mu) for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 20,0,TMath::Pi());
	TH1D *histTroubleDEtaMuMu    = new TH1D("hTroubleDEtaMuMu","#Delta #eta (#mu-#mu) for trks #Delta R(#mu-trk) = 0.7 - 1.2 ; Value; N_{events}", 20,-5,5);
	TH1D *histTroubleMu1Pt       = new TH1D("hTroubleMu1Pt","Mu1 p_{T} for trks #Delta R(#mu1-trk) = 0.7 - 1.2 ; Value; N_{events}", 25,0,35);
	TH1D *histTroubleMu2Pt       = new TH1D("hTroubleMu2Pt","Mu2 p_{T} for trks #Delta R(#mu2-trk) = 0.7 - 1.2 ; Value; N_{events}", 20,0,25);
	TH2D *histTroubleEtaVsPhi1   = new TH2D("hTroubleEtaVsPhi1","dPhi vs dEta of tracks (>2.5 GeV) vs muon 1 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());
	TH2D *histTroubleEtaVsPhi2   = new TH2D("hTroubleEtaVsPhi2","dPhi vs dEta of tracks (>2.5 GeV) vs muon 2 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// Plots for testing invariant mass correlation
	double massBins[6]           = {0,1,2,3,4,10};
	TH1D *histM1                 = new TH1D("hM1", "Inv. Mass of 1st system, full selection; m(#mu_{1}-tk) [GeV]; N_{events}", 10,0,10);
	TH1D *histM2                 = new TH1D("hM2", "Inv. Mass of 2st system, full selection; m(#mu_{2}-tk) [GeV]; N_{events}", 10,0,10);

	// MC truth - use actual mu-tk pairs from tau
	TH1D *histM1_truth_0to1      = new TH1D("hM1_truth_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_1to2      = new TH1D("hM1_truth_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_2to3      = new TH1D("hM1_truth_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_truth_3toInf    = new TH1D("hM1_truth_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	
	// actual dist using signal selection
	TH1D *histM1_0to1            = new TH1D("hM1_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_1to2            = new TH1D("hM1_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_2to3            = new TH1D("hM1_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_3toInf          = new TH1D("hM1_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);

	// actual dist using sideband selection
	TH1D *histM1_side_0to1       = new TH1D("hM1_side_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_side_1to2       = new TH1D("hM1_side_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_side_2to3       = new TH1D("hM1_side_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D *histM1_side_3toInf     = new TH1D("hM1_side_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);

	int nMu(0);
	int n2p5(0), n2p5OS(0); // count # muons with 1+ tracks with pT > 2.5 for SS+OS, and OS
	int nMuPass(0);

	// Loop over all events
	Long64_t numberOfEntries = treeReader->GetEntries();
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
		// Track *candTk(0);

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

			histHPt->Fill(((a1->P4())+(a2->P4())).Pt());

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
			histDRMuMu->Fill(mu1Mom.DeltaR(mu2Mom));
			histDEtaVsDPhiMuMu->Fill(fabs(origMu1->Eta-origMu2->Eta),fabs(mu1Mom.DeltaPhi(mu2Mom)));

			histMu1PtSel->Fill(mu1PT);
			histMu2PtSel->Fill(mu2PT);

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
			bool atLeastTk2p5 = false; // to monitor if theres a tk with pT > 2.5
			bool atLeastTk2p5OS = false; // same but for OS tk-muon
			for(int a = 0; a < branchTracks->GetEntries(); a++){
				candTk = (Track*) branchTracks->At(a);

				if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
					&& (candTk->PT != mu2->PT)
					&& (candTk->PT > 1.)
					&& (fabs(candTk->Z) < 1.) //dz < 1mm
					&& ((pow(candTk->X,2)+pow(candTk->Y,2)) < 1.) //dxy < 1mm
					&& (fabs(candTk->Eta)<3)
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

						histNTracks1->Fill(dR1);
						histNTracks2->Fill(dR2);
						atLeastTk2p5 = true;

						histTroubleEtaVsPhi1->Fill(fabs(candTk->Eta - mu1Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu1Mom)));
						histTroubleEtaVsPhi2->Fill(fabs(candTk->Eta - mu2Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu2Mom)));

						if (dR1 < 0.5){
							tk1_2p5.push_back(candTk);
						}
						if (dR2 < 0.5){
							tk2_2p5.push_back(candTk);
						}
						if ((candTk->Charge) * (mu1->Charge) < 0) {
							histNTracks1OS->Fill(dR1);
							histNTracks2OS->Fill(dR2);
							atLeastTk2p5OS = true;

							if (dR1 < 0.5){
								tk1_2p5_OS.push_back(candTk);
							}
							if (dR2 < 0.5){
								tk2_2p5_OS.push_back(candTk);
							}
						}					
					}  else {
						if (dR1 < 0.5){
							tk1_1to2p5.push_back(candTk);
						}
						if (dR2 < 0.5){
							tk2_1to2p5.push_back(candTk);
						}
					}
				} // End of track selection
			} // End of track loop

			// Count # muons that contribute to track distribution plots
			if (atLeastTk2p5) n2p5++;
			if (atLeastTk2p5OS) n2p5OS++;

			// SIGNAL SELECTION
			if (tk1_1.size() == 1 && tk2_1.size() == 1 
			&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {
				
				nMuPass++;
				histTrack1Pt->Fill(tk1_2p5_OS[0]->PT);
				histTrack2Pt->Fill(tk2_2p5_OS[0]->PT);

				TLorentzVector track1Mom=tk1_2p5_OS[0]->P4();
				TLorentzVector track2Mom=tk2_2p5_OS[0]->P4();

				// combined mu+tk system
				TLorentzVector sys1 = mu1Mom+track1Mom;
				TLorentzVector sys2 = mu2Mom+track2Mom;

				// plot the pT, dR, dEta, DPhi of two systems
				histSys1Pt->Fill(sys1.Pt());
				histSys2Pt->Fill(sys2.Pt());
				histDRSys->Fill(sys1.DeltaR(sys2));
				histDEtaVsDPhiSys->Fill(fabs(sys1.Eta() - sys2.Eta()),fabs(sys1.DeltaPhi(sys2)));

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

		// // SIDEBAND REGION
		// // one muon has 1 tk > 2.5 (OS), other has 2 or 3 ( 1 tk => "mu1", 2/3 tk => "mu2")
		// // so NOT pT ordered (although you could do that by eliminating the else if ... bit)
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

	// Clone and rescale some hists
	TH1D* histNTracksAbs1   = (TH1D*)histNTracks1->Clone("hNTracksAbs1");
	TH1D* histNTracksAbs2   = (TH1D*)histNTracks2->Clone("hNTracksAbs2");
	TH1D* histNTracksAbs1OS = (TH1D*)histNTracks1OS->Clone("hNTracksAbs1OS");
	TH1D* histNTracksAbs2OS = (TH1D*)histNTracks2OS->Clone("hNTracksAbs2OS");

	// Rescale some hists
	// Abs # of tracks per muon
	histNTracksAbs1->Scale(1./n2p5);
	histNTracksAbs1OS->Scale(1./n2p5OS);
	histNTracksAbs2->Scale(1./n2p5);
	histNTracksAbs2OS->Scale(1./n2p5OS);

	histNTracksAbs1->SetYTitle("Ave. N_{trk} per #mu_{1}");
	histNTracksAbs2->SetYTitle("Ave. N_{trk} per #mu_{2}");
	histNTracksAbs1OS->SetYTitle("Ave. OS N_{trk} per #mu_{1}");
	histNTracksAbs2OS->SetYTitle("Ave. OS N_{trk} per #mu_{2}");

	// AU scaling
	histNTracks1->Scale(1./histNTracks1->Integral());
	histNTracks2->Scale(1./histNTracks2->Integral());
	histNTracks1OS->Scale(1./histNTracks1OS->Integral());
	histNTracks2OS->Scale(1./histNTracks2OS->Integral());

	// Cumulative plots
	histNTracksCum1   = (TH1D*)histNTracksAbs1->Clone("hNTracksCum1");
	histNTracksCum2   = (TH1D*)histNTracksAbs2->Clone("hNTracksCum2");
	histNTracksCum1OS = (TH1D*)histNTracksAbs1OS->Clone("hNTracksCum1OS");
	histNTracksCum2OS = (TH1D*)histNTracksAbs2OS->Clone("hNTracksCum2OS");

	for (int i = 1; i <= histNTracks1->GetNbinsX(); i++){
		histNTracksCum1->SetBinContent(i,histNTracksCum1->GetBinContent(i-1) + histNTracksAbs1->GetBinContent(i));
		histNTracksCum2->SetBinContent(i,histNTracksCum2->GetBinContent(i-1) + histNTracksAbs2->GetBinContent(i));
		histNTracksCum1OS->SetBinContent(i,histNTracksCum1OS->GetBinContent(i-1) + histNTracksAbs1OS->GetBinContent(i));
		histNTracksCum2OS->SetBinContent(i,histNTracksCum2OS->GetBinContent(i-1) + histNTracksAbs2OS->GetBinContent(i));
	}

	cout << "# muons with 1+ track with pT >2.5 GeV: " << n2p5 << endl;
	cout << "# muons with 1+ OS track with pT > 2.5GeV: " << n2p5OS << endl;
	cout << "nMuPass: " << nMuPass << endl;

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

	// drawHistAndSave(histNMu, "HISTE", "NMu", directory, app);
	// drawHistAndSave(histNMu1, "HISTE", "NMu1", directory, app);
	// drawHistAndSave(histNMu2, "HISTE", "NMu2", directory, app);

	drawHistAndSave(histMu1Pt, "HISTE", "Mu1Pt", directory, app);
	drawHistAndSave(histMu2Pt, "HISTE", "Mu2Pt", directory, app);
	drawHistAndSave(histTrack1Pt, "HISTE", "Track1Pt", directory, app);
	drawHistAndSave(histTrack2Pt, "HISTE", "Track2Pt", directory, app);

	drawHistAndSave(histSys1Pt, "HISTE", "Sys1Pt", directory, app);
	drawHistAndSave(histSys2Pt, "HISTE", "Sys2Pt", directory, app);
	drawHistAndSave(histDRSys, "HISTE", "DRSys", directory, app);
	drawHistAndSave(histDEtaVsDPhiSys, "COLZ", "DEtaVsDPhiSys", directory, app);

	drawHistAndSave(histSys1PtTruth, "HISTE", "Sys1PtTruth", directory, app);
	drawHistAndSave(histSys2PtTruth, "HISTE", "Sys2PtTruth", directory, app);
	drawHistAndSave(histDRSysTruth, "HISTE", "DRSysTruth", directory, app);
	drawHistAndSave(histDEtaVsDPhiSysTruth, "COLZ", "DEtaVsDPhiSysTruth", directory, app);

	drawHistAndSave(histNuPt, "HISTE", "NuPt", directory, app);
	
	drawHistAndSave(histMu1PtSel, "HISTE", "Mu1PtSel", directory, app);
	drawHistAndSave(histMu2PtSel, "HISTE", "Mu2PtSel", directory, app);

	// Track distributions around muons
	drawHistAndSave(histNTracks1, "HISTE", "NTracks1_NS", directory, app);
	drawHistAndSave(histNTracks2, "HISTE", "NTracks2_NS", directory, app);
	drawHistAndSave(histNTracks1OS, "HISTE", "NTracks1_OS", directory, app);
	drawHistAndSave(histNTracks2OS, "HISTE", "NTracks2_OS", directory, app);

	drawHistAndSave(histNTracksAbs1, "HISTE", "NTracksAbs1_NS", directory, app);
	drawHistAndSave(histNTracksAbs2, "HISTE", "NTracksAbs2_NS", directory, app);
	drawHistAndSave(histNTracksAbs1OS, "HISTE", "NTracksAbs1_OS", directory, app);
	drawHistAndSave(histNTracksAbs2OS, "HISTE", "NTracksAbs2_OS", directory, app);

	drawHistAndSave(histNTracksCum1, "HISTE", "NTracksCum1_NS", directory, app);
	drawHistAndSave(histNTracksCum2, "HISTE", "NTracksCum2_NS", directory, app);
	drawHistAndSave(histNTracksCum1OS, "HISTE", "NTracksCum1_OS", directory, app);
	drawHistAndSave(histNTracksCum2OS, "HISTE", "NTracksCum2_OS", directory, app);

	drawHistAndSave(histDRMuMu, "HISTE", "DRMuMu", directory, app);
	drawHistAndSave(histDEtaVsDPhiMuMu, "COLZ", "DEtaVsDPhiMuMu", directory, app);

	drawHistAndSave(histNTk, "HISTE", "NTk", directory, app);
	drawHistAndSave(histNTk1, "HISTE", "NTk1", directory, app);
	drawHistAndSave(histNTk25, "HISTE", "NTk25", directory, app);

	if (doSignal){
		drawHistAndSave(histDRa1, "HISTE", "DRa1", directory, app);
		drawHistAndSave(histDRa2, "HISTE", "DRa2", directory, app);
		drawHistAndSave(histPID, "HISTE", "PID", directory, app);
		drawHistAndSave(histHPt, "HISTE", "HPt", directory, app);
	}

	drawHistAndSave(histTroublePt, "HISTE", "TroubleTkPt", directory, app);
	drawHistAndSave(histTroublePID, "HISTE", "TroubleTkPID", directory, app);
	drawHistAndSave(histTroubleEta, "HISTE", "TroubleTkEta", directory, app);
	drawHistAndSave(histTroublePhi, "HISTE", "TroubleTkPhi", directory, app);
	drawHistAndSave(histTroubleMatch, "HISTE", "TroubleMatch", directory, app);
	drawHistAndSave(histTroubleDRMuMu, "HISTE", "TroubleDRMuMu", directory, app);

	drawHistAndSave(histTroubleDPhiMuMu, "HISTE", "TroubleDPhiMuMu", directory, app);
	drawHistAndSave(histTroubleDEtaMuMu, "HISTE", "TroubleDEtaMuMu", directory, app);
	drawHistAndSave(histTroubleMu1Pt, "HISTE", "TroubleMu1Pt", directory, app);
	drawHistAndSave(histTroubleMu2Pt, "HISTE", "TroubleMu2Pt", directory, app);

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
	c.SaveAs((directory+"/TroubleEtaVsPhi1_"+delph+"_"+app+".pdf").c_str());
	histTroubleEtaVsPhi2->Draw("COLZ");
	problemRing1.Draw("only");
	problemRing2.Draw("only");
	c.SaveAs((directory+"/TroubleEtaVsPhi2_"+delph+"_"+app+".pdf").c_str());

	drawHistAndSave(histRand, "HISTE", "RandTest", directory, app);

	// Mass plots
	drawHistAndSave(histM1, "HISTE", "M1", directory, app);
	drawHistAndSave(histM2, "HISTE", "M2", directory, app);
	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - signal region;m(#mu_{1}-tk) [GeV]; A.U.", histM1_0to1, histM1_1to2, histM1_2to3, histM1_3toInf, "M1_M2", directory, app);
	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - sideband region;m(#mu_{1}-tk) [GeV]; A.U.", histM1_side_0to1, histM1_side_1to2, histM1_side_2to3, histM1_side_3toInf, "M1_M2_side", directory, app);
	if(doSignal){
		drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - MC truth;m(#mu_{1}-tk) [GeV]; A.U.", histM1_truth_0to1, histM1_truth_1to2, histM1_truth_2to3, histM1_truth_3toInf, "M1_M2_truth", directory, app);
	}

	TFile* outFile = TFile::Open((directory+"/output_"+delph+"_"+app+".root").c_str(),"RECREATE");

	histMu1Pt->Write("",TObject::kOverwrite);
	histMu2Pt->Write("",TObject::kOverwrite);
	histTrack1Pt->Write("",TObject::kOverwrite);
	histTrack2Pt->Write("",TObject::kOverwrite);
	histSys1Pt->Write("",TObject::kOverwrite);
	histSys2Pt->Write("",TObject::kOverwrite);
	histDRSys->Write("",TObject::kOverwrite);
	histDEtaVsDPhiSys->Write("",TObject::kOverwrite);
	histSys1PtTruth->Write("",TObject::kOverwrite);
	histSys2PtTruth->Write("",TObject::kOverwrite);
	histDRSysTruth->Write("",TObject::kOverwrite);
	histDEtaVsDPhiSysTruth->Write("",TObject::kOverwrite);
	histNuPt->Write("",TObject::kOverwrite);
	histMu1PtSel->Write("",TObject::kOverwrite);
	histMu2PtSel->Write("",TObject::kOverwrite);
	histNTracks1->Write("",TObject::kOverwrite);
	histNTracks2->Write("",TObject::kOverwrite);
	histNTracks1OS->Write("",TObject::kOverwrite);
	histNTracks2OS->Write("",TObject::kOverwrite);
	histNTracksAbs1->Write("",TObject::kOverwrite);
	histNTracksAbs2->Write("",TObject::kOverwrite);
	histNTracksAbs1OS->Write("",TObject::kOverwrite);
	histNTracksAbs2OS->Write("",TObject::kOverwrite);
	histNTracksCum1->Write("",TObject::kOverwrite);
	histNTracksCum2->Write("",TObject::kOverwrite);
	histNTracksCum1OS->Write("",TObject::kOverwrite);
	histNTracksCum2OS->Write("",TObject::kOverwrite);

	histDRMuMu->Write("",TObject::kOverwrite);
	histDEtaVsDPhiMuMu->Write("",TObject::kOverwrite);
	histNTk->Write("",TObject::kOverwrite);
	histNTk1->Write("",TObject::kOverwrite);
	histNTk25->Write("",TObject::kOverwrite);

	if (doSignal){
		histDRa1->Write("",TObject::kOverwrite);
		histDRa2->Write("",TObject::kOverwrite);
		histPID->Write("",TObject::kOverwrite);
	}

	histTroublePt->Write("",TObject::kOverwrite);
	histTroublePID->Write("",TObject::kOverwrite);
	histTroubleEta->Write("",TObject::kOverwrite);
	histTroublePhi->Write("",TObject::kOverwrite);
	histTroubleMatch->Write("",TObject::kOverwrite);
	histTroubleDRMuMu->Write("",TObject::kOverwrite);
	histTroubleDPhiMuMu->Write("",TObject::kOverwrite);
	histTroubleDEtaMuMu->Write("",TObject::kOverwrite);
	histTroubleMu1Pt->Write("",TObject::kOverwrite);
	histTroubleMu2Pt->Write("",TObject::kOverwrite);
	// Mass plots
	histM1->Write("",TObject::kOverwrite);
	histM2->Write("",TObject::kOverwrite);

	histM1_truth_0to1->Write("",TObject::kOverwrite);
	histM1_truth_1to2->Write("",TObject::kOverwrite);
	histM1_truth_2to3->Write("",TObject::kOverwrite);
	histM1_truth_3toInf->Write("",TObject::kOverwrite);
	histM1_0to1->Write("",TObject::kOverwrite);
	histM1_1to2->Write("",TObject::kOverwrite);
	histM1_2to3->Write("",TObject::kOverwrite);
	histM1_3toInf->Write("",TObject::kOverwrite);
	histM1_side_0to1->Write("",TObject::kOverwrite);
	histM1_side_1to2->Write("",TObject::kOverwrite);
	histM1_side_2to3->Write("",TObject::kOverwrite);
	histM1_side_3toInf->Write("",TObject::kOverwrite);

	outFile->Close();

}
