#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include "commonFunctions.h"
#include "cuts.h"

using std::cout;
using std::endl;

/**
 * Main analysis script to make plots, except for those pertianing to mass distributions or correlation coefficiants (see massPlots.C) 
 */

/**
 * makes clone fo hist, appends "suffix" to hist name, 
 * write to currently open file
 */
template <typename T>
void makeCopySave(T* h, std::string suffix="_unnormalised") {
	T* h_clone = (T*) h->Clone((h->GetName()+suffix).c_str());
	h_clone->Write("", TObject::kOverwrite);
}

void mainAnalysis(int argc, char* argv[])
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");
	gStyle->SetHistLineWidth(2);

	ProgramOpts pOpts(argc, argv);

	// MCsource source     = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	bool doSignal       = pOpts.getSignal(); // for signal or not
	// bool doMu           = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = pOpts.getMuOrdering(); // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT          = pOpts.getHLT(); // whether to use MC that has HLT cuts already applied or not.
	double deltaR       = pOpts.getdR(); // dR(mu-mu) value to use

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
	// TClonesArray *branchGenMuons = treeReader->UseBranch("OnlyGenMuons"); // GenParticle object
	TClonesArray *branchGenMuons = treeReader->UseBranch("GenMuon"); // Track object
	// TClonesArray *branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray *branchAll      = treeReader->UseBranch("AllParticle");

	//////////////////////
	// Book histograms  //
	//////////////////////
	
	// Track distributions (pT > 2.5)
	TH1D *histNTracks1           = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks2           = new TH1D("hNTracks2" ,"Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);
	TH1D *histNTracks1OS         = new TH1D("hNTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks2OS         = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);

	TH1D *histNTracksCum1        = new TH1D("hNTracksCum1" ,"Cumulative Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNTracksCum2        = new TH1D("hNTracksCum2" ,"Cumulative Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);
	TH1D *histNTracksCum1OS      = new TH1D("hNTracksCum1OS" ,"Cumulative Number of tracks about mu1, OS,p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNTracksCum2OS      = new TH1D("hNTracksCum2OS" ,"Cumulative Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);
	
	// Track distributions (pT > 1)
	TH1D *histNTracksAll1           = new TH1D("hNTracksAll1" ,"Number of tracks about mu1, p_{T}(trk)>1 GeV;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracksAll2           = new TH1D("hNTracksAll2" ,"Number of tracks about mu2, p_{T}(trk)>1 GeV;#Delta R (#mu_{2}-track); A.U.", 50,0,5);


	// Soft track distributions, for events where you already have 1 track with pT > 2.5
	TH1D *histNSoftTracks1       = new TH1D("hNSoftTracks1" ,"Number of tracks about mu1, p_{T}(trk) = (1-2.5)GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNSoftTracks2       = new TH1D("hNSoftTracks2" ,"Number of tracks about mu2, p_{T}(trk) = (1-2.5)GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{2}-track); A.U.", 50,0,5);
	TH1D *histNSoftTracks1OS     = new TH1D("hNSoftTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk) = (1-2.5)GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNSoftTracks2OS     = new TH1D("hNSoftTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk) = (1-2.5)GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{2}-track); A.U.", 50,0,5);

	TH1D *histNSoftTracksCum1    = new TH1D("hNSoftTracksCum1" ,"Cumulative Number of tracks about mu1, p_{T}(trk) = (1-2.5) GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNSoftTracksCum2    = new TH1D("hNSoftTracksCum2" ,"Cumulative Number of tracks about mu2, p_{T}(trk) = (1-2.5) GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);
	TH1D *histNSoftTracksCum1OS  = new TH1D("hNSoftTracksCum1OS" ,"Cumulative Number of tracks about mu1, OS,p_{T}(trk) = (1-2.5) GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNSoftTracksCum2OS  = new TH1D("hNSoftTracksCum2OS" ,"Cumulative Number of tracks about mu2, OS, p_{T}(trk) = (1-2.5) GeV, muon selection, 1 OS trk with pT > 2.5 GeV;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);

	// Particle kinematics
	TH1D *histMu1Pt              = new TH1D("hMu1Pt", "#mu_{1} p_{T}, no selection ;#mu_{1} p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histMu2Pt              = new TH1D("hMu2Pt", "#mu_{2} p_{T}, no selection;#mu_{2} p_{T} [GeV]; N_{events}", 50,0,50.);

	TH1D *histTrack1Pt           = new TH1D("hTrack1Pt","Track 1 p_{T}, signal selection; track 1 p_{T} [GeV]; N_{events}", 25,0,25.);
	TH1D *histTrack2Pt           = new TH1D("hTrack2Pt","Track 2 p_{T}, signal selection; track 2 p_{T} [GeV]; N_{events}", 25,0,25.);
	
	// For combined mu+tk systems
	TH1D *histSys1Pt             = new TH1D("hSys1Pt", "System 1 p_{T}, signal selection ;System 1 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histSys2Pt             = new TH1D("hSys2Pt", "System 2 p_{T}, signal selection;System 2 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histDRSys              = new TH1D("hDRSys", "#Delta R(Sys_{1}-Sys_{2}), signal selection;#Delta R(Sys_{1}-Sys_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiSys      = new TH2D("hDEtaVsDPhiSys","dPhi vs dEta for system 2 wrt system 1 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// For combined mu+tk systems - MC Truth
	TH1D *histSys1PtTruth        = new TH1D("hSys1PtTruth", "System 1 p_{T}, signal selection, MC truth;System 1 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histSys2PtTruth        = new TH1D("hSys2PtTruth", "System 2 p_{T}, signal selection, MC truth;System 2 p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histDRSysTruth         = new TH1D("hDRSysTruth", "#Delta R(Sys_{1}-Sys_{2}), signal selection, MC truth;#Delta R(Sys_{1}-Sys_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiSysTruth = new TH2D("hDEtaVsDPhiSysTruth","dPhi vs dEta for system 2 wrt system 1, MC Truth ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	
	TH1D *histMu1PtSel           = new TH1D("hMu1PtSel", "#mu_{1} p_{T}, muon selection, no tk selection;#mu_{1} p_{T} [GeV]; N_{events}", 50,0,50.);
	TH1D *histMu2PtSel           = new TH1D("hMu2PtSel", "#mu_{2} p_{T}, muon selection, no tk selection;#mu_{2} p_{T} [GeV]; N_{events}", 50,0,50.);

	TH1D *histDRMuMu             = new TH1D("hDRMuMu", "#Delta R(#mu-#mu), muon selection;#Delta R(#mu_{1}-#mu_{2}); N_{events}", 30,0,5);
	TH2D *histDEtaVsDPhiMuMu     = new TH2D("hDEtaVsDPhiMuMu","dPhi vs dEta of selection muons ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// Signal-specific ones
	TH1D *histDRa1               = new TH1D("hDRa1","#Delta R(#tau-#tau) 1st a_{1}, no muon selection;#Delta R(#tau-#tau); N_{events}", 10,0,1.);
	TH1D *histDRa2               = new TH1D("hDRa2","#Delta R(#tau-#tau) 2nd a_{1}, no muon selection;#Delta R(#tau-#tau); N_{events}", 10,0,1.);

	TH1D *histHPt                = new TH1D("hHPt","h_{1,2} p_{T}, no selection ;h_{1,2} p_{T} [GeV]; N_{events}",25,0,50); // Isn't included in MC file from Calchep :(
	TH1D *histNuPt               = new TH1D("hNuPt", "#nu p_{T}, no selection ;#nu p_{T} [GeV]; N_{events}", 50,0,50.);

	TH1D *histPID                = new TH1D("hPID","PID of tau 1-prong; PID; N_{events}", 350,0,350);

	// TH1D *histRand            = new TH1D("hRand","Testing rand(); Value; N_{events}", 100,0,1);
	
	TH2D *histTkEtaVsPhi1        = new TH2D("hTkEtaVsPhi1","dPhi vs dEta of tracks (>2.5 GeV) vs muon 1 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());
	TH2D *histTkEtaVsPhi2        = new TH2D("hTkEtaVsPhi2","dPhi vs dEta of tracks (>2.5 GeV) vs muon 2 ; #Delta #eta; #Delta #phi", 30,0,3, 20, 0, TMath::Pi());

	// int nMu(0);
	int n2p5(0), n2p5OS(0); // count # muons with 1+ tracks with pT > 2.5 for SS+OS, and OS
	int nOnly2p5(0), nOnly2p5OS(0); // count # muons with 1 tracks with  pT > 2.5 for SS+OS, and OS
	int n1(0);
	int nMuPass(0);

	//////////////////////
	// Loop over events //
	//////////////////////
	Long64_t numberOfEntries = getNumberEvents(treeReader, &pOpts);
	cout << "Running over " << numberOfEntries << " events" << endl;

	bool stop = false; // used to stop the loop, for debugging/testing

	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		// Do at gen particle level.
		// histNMu->Fill(branchGenMuons->GetEntries());

		if (branchGenMuons->GetEntries() < 2) continue; // skip if <2 muons!

		//////////////////////////////////////////////////
		// Get the hard interaction particles for signal MC //
		// No selection cuts applied (only >=2 muons)       //
		//////////////////////////////////////////////////////
		
		GenParticle *charged1a(nullptr);
		GenParticle *charged1b(nullptr);
		GenParticle *charged2a(nullptr);
		GenParticle *charged2b(nullptr);

		if (doSignal) {
			GenParticle *a1(nullptr), *a2(nullptr), *cand(nullptr);
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
						// cout << "found first a1 at " << j << " check PID " << a1->PID << endl;
					} else {
						// cout << "found second a1 at " << j << endl;
						a2=cand;
					}
				}
			}

			histHPt->Fill(((a1->P4())+(a2->P4())).Pt());

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
				GenParticle* muTruth1(nullptr);
				GenParticle* trackTruth1(nullptr);
				GenParticle* muTruth2(nullptr);
				GenParticle* trackTruth2(nullptr);

				// Assign charged products to be mu or track
				bool truth1HasMu = assignMuonAndTrack(muTruth1, trackTruth1, *charged1a, *charged1b);				
				bool truth2HasMu = assignMuonAndTrack(muTruth2, trackTruth2, *charged2a, *charged2b);

				// NOTE: muons are NOT pT ordered

				if (!truth1HasMu || !truth2HasMu) {
					cout << "Problem, no truth mu for 1 and/or 2!" << endl;
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

				// Dirty!
				if (!muTruth1) delete muTruth1;
				if (!trackTruth1) delete trackTruth1;
				if (!muTruth2) delete muTruth2;
				if (!trackTruth2) delete trackTruth2;
			} // end if(charged1a...) 
			else {
				throw runtime_error("Not all prongs found!");
			}
		} // end if(doSignal)

		/////////////////////////////////////////////////////////////////////////
		// Now, get the two highest pT muons in the event that pass selection, // 
		// store pointers to the Track particles and 4-momenta                 //
		// (Use tracks for muons as store more info about position)
		/////////////////////////////////////////////////////////////////////////
		
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


		// Make pairs, see if they pass all cuts (SS, eta, deltaR, dZ, d0)
		// If they do, store in mu1 and mu2 (mu1 has higher pT)
		std::pair<Track*, Track*> p = testMuons(muons17toInf, muons10to17, &checkMuons, deltaR);
		Track* mu1 = p.first;
		Track* mu2 = p.second;

		if (!(p.first && p.second)) continue;

		// Now randomly swap mu1 - mu2
		Track *origMu1(nullptr), *origMu2(nullptr);
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

		// Plot some muon quantities
		histMu1Pt->Fill(mu1->PT);
		histMu2Pt->Fill(mu2->PT);
		histDRMuMu->Fill(mu1Mom.DeltaR(mu2Mom));
		histDEtaVsDPhiMuMu->Fill(fabs(origMu1->Eta-origMu2->Eta),fabs(mu1Mom.DeltaPhi(mu2Mom)));
		histMu1PtSel->Fill(mu1->PT);
		histMu2PtSel->Fill(mu2->PT);

		/////////////////////////////////
		// Look at tracks around muons //
		/////////////////////////////////

		// Vectors of tracks with pT > 1, within dR < 0.5 of respective muons + other cuts
		// so tk1 is the track nearest to muon1, (may or may not be highest pT, depends if random swapping is on)
		std::vector<Track*> tk1_1;
		std::vector<Track*> tk2_1;

		// same but all dR
		std::vector<Track*> tk1_1_alldR;
		std::vector<Track*> tk2_1_alldR;

		// same but with pT >2.5
		std::vector<Track*> tk1_2p5;
		std::vector<Track*> tk2_2p5;

		// same but with pT > 2.5, OS to muon
		std::vector<Track*> tk1_2p5_OS;
		std::vector<Track*> tk2_2p5_OS;
		
		// same but with 1 < pT < 2.5 (for sideband)
		std::vector<Track*> tk1_1to2p5;
		std::vector<Track*> tk2_1to2p5;

		// same but with 1 < pT < 2.5, no dR restriction 
		std::vector<Track*> tk1_1to2p5_alldR;
		std::vector<Track*> tk2_1to2p5_alldR;
		
		Track *candTk(nullptr);
		bool atLeastTk1 = false;
		bool atLeastTk2p5 = false; // to monitor if theres a tk with pT > 2.5
		bool atLeastTk2p5OS = false; // same but for OS tk-muon
		for(int a = 0; a < branchTracks->GetEntries(); a++){
			candTk = (Track*) branchTracks->At(a);

			if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
				&& (candTk->PT != mu2->PT)
				&& checkTrackLoose(candTk)
			){
				// Store track in suitable vector
				double dR1 = (candTk->P4()).DeltaR(mu1Mom);
				double dR2 = (candTk->P4()).DeltaR(mu2Mom);

				fillTrackVectors(candTk, mu1, mu2, &tk1_1, &tk2_1);
				fillTrackVectors(candTk, mu1, mu2, &tk1_1_alldR, &tk2_1_alldR, 5);

				if (checkTrackTight(candTk)){

					histNTracks1->Fill(dR1);
					histNTracks2->Fill(dR2);
					atLeastTk2p5 = true;

					histTkEtaVsPhi1->Fill(fabs(candTk->Eta - mu1Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu1Mom)));
					histTkEtaVsPhi2->Fill(fabs(candTk->Eta - mu2Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu2Mom)));

					fillTrackVectors(candTk, mu1, mu2, &tk1_2p5, &tk2_2p5);

					if (checkTkMuOS(candTk, mu1)) {
						histNTracks1OS->Fill(dR1);
						histNTracks2OS->Fill(dR2);
						atLeastTk2p5OS = true;

						fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS, &tk2_2p5_OS);

					}					
				}  else {
					tk1_1to2p5_alldR.push_back(candTk);
					tk2_1to2p5_alldR.push_back(candTk);
					fillTrackVectors(candTk, mu1, mu2, &tk1_1to2p5, &tk2_1to2p5);

				}
			} // End of track selection
		} // End of track loop

		// Count # muons that contribute to track distribution plots
		// if (atLeastTk1) n1++;
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
				
				// For only if at least one muon has a soft track
				if(!(tk1_1to2p5.size() == 0 && tk2_1to2p5.size() == 0)){ 
				}

		} // end of sideband 1to2p5 

		// Slightly different region - for additional track investigations
		// For soft track distributions
		// ATM it uses signal region. Not put in the signal region bit above, as subject to future modification
		if (tk1_1.size() >= 1 && tk2_1.size() >= 1 
			&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1){
			
			nOnly2p5OS++;
			if (tk1_1to2p5_alldR.size() > 0){
				for( auto softTk : tk1_1to2p5_alldR ){
					histNSoftTracks1->Fill((softTk->P4()).DeltaR(mu1Mom));
					if (softTk->Charge * mu1->Charge < 0)
						histNSoftTracks1OS->Fill((softTk->P4()).DeltaR(mu1Mom));
				}
			}
			if (tk2_1to2p5_alldR.size() > 0){
				for( auto softTk : tk2_1to2p5_alldR ){
					histNSoftTracks2->Fill((softTk->P4()).DeltaR(mu2Mom));
					if (softTk->Charge * mu2->Charge < 0)
						histNSoftTracks2OS->Fill((softTk->P4()).DeltaR(mu2Mom));
				}
			}
		}

		// another sideband to plot all tracks with pT > 1
		if (tk1_1.size() >= 1 && tk2_1.size() >= 1 	
			&& tk1_2p5_OS.size() >= 1 && tk2_2p5_OS.size() >= 1){
			
			n1++;
			if (tk1_1_alldR.size() > 0){
				for (auto tk: tk1_1_alldR) {
					histNTracksAll1->Fill(tk->P4().DeltaR(mu1Mom));
				}
			}
			if (tk2_1_alldR.size() > 0){
				for (auto tk: tk2_1_alldR) {
					histNTracksAll2->Fill(tk->P4().DeltaR(mu1Mom));
				}
			}
		}

	} // end of event loop

	// Clone and rescale some hists
	TH1D* histNTracksAbs1   = (TH1D*)histNTracks1->Clone("hNTracksAbs1");
	TH1D* histNTracksAbs2   = (TH1D*)histNTracks2->Clone("hNTracksAbs2");
	TH1D* histNTracksAbs1OS = (TH1D*)histNTracks1OS->Clone("hNTracksAbs1OS");
	TH1D* histNTracksAbs2OS = (TH1D*)histNTracks2OS->Clone("hNTracksAbs2OS");

	TH1D* histNSoftTracksAbs1   = (TH1D*)histNSoftTracks1->Clone("hNSoftTracksAbs1");
	TH1D* histNSoftTracksAbs2   = (TH1D*)histNSoftTracks2->Clone("hNSoftTracksAbs2");
	TH1D* histNSoftTracksAbs1OS = (TH1D*)histNSoftTracks1OS->Clone("hNSoftTracksAbs1OS");
	TH1D* histNSoftTracksAbs2OS = (TH1D*)histNSoftTracks2OS->Clone("hNSoftTracksAbs2OS");

	TH1D* histNTracksAllAbs1   = (TH1D*)histNTracksAll1->Clone("hNTracksAllAbs1");
	TH1D* histNTracksAllAbs2   = (TH1D*)histNTracksAll2->Clone("hNTracksAllAbs2");

	// Rescale some hists
	// Abs # of tracks per muon
	histNTracksAbs1->Scale(1./n2p5);
	histNTracksAbs1OS->Scale(1./n2p5OS);
	histNTracksAbs2->Scale(1./n2p5);
	histNTracksAbs2OS->Scale(1./n2p5OS);
	
	histNTracksAllAbs1->Scale(1./n1);
	histNTracksAllAbs2->Scale(1./n1);
	histNTracksAllAbs1->SetYTitle("Ave. N_{trk} per #mu_{1}");
	histNTracksAllAbs2->SetYTitle("Ave. N_{trk} per #mu_{2}");

	histNTracksAbs1->SetYTitle("Ave. N_{trk} per #mu_{1}");
	histNTracksAbs2->SetYTitle("Ave. N_{trk} per #mu_{2}");
	histNTracksAbs1OS->SetYTitle("Ave. OS N_{trk} per #mu_{1}");
	histNTracksAbs2OS->SetYTitle("Ave. OS N_{trk} per #mu_{2}");

	histNSoftTracksAbs1->Scale(1./nOnly2p5OS);
	histNSoftTracksAbs1OS->Scale(1./nOnly2p5OS);
	histNSoftTracksAbs2->Scale(1./nOnly2p5OS);
	histNSoftTracksAbs2OS->Scale(1./nOnly2p5OS);

	histNSoftTracksAbs1->SetYTitle("Ave. N_{trk} per #mu_{1}");
	histNSoftTracksAbs2->SetYTitle("Ave. N_{trk} per #mu_{2}");
	histNSoftTracksAbs1OS->SetYTitle("Ave. OS N_{trk} per #mu_{1}");
	histNSoftTracksAbs2OS->SetYTitle("Ave. OS N_{trk} per #mu_{2}");


	////////////////////////
	// Draw hists to PDF //
	////////////////////////

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

	app += "_dR";
	app += boost::lexical_cast<std::string>(deltaR);

	// Get directory that input file was in - put plots in there
	std::string directory = getDirectory(chain.GetFile());
	// Get Delphes file config used - last part of directory name
	std::string delph = getDelph(directory);

	TFile* outFile = TFile::Open((directory+"/output_main_"+delph+"_"+app+".root").c_str(),"UPDATE");
	cout << "Writing to " << outFile->GetName() << endl;

	// AU scaling
	// save unnormalised version to file
	makeCopySave(histNTracks1);
	normaliseHist(histNTracks1);
	makeCopySave(histNTracks2);
	normaliseHist(histNTracks2);
	makeCopySave(histNTracks1OS);
	normaliseHist(histNTracks1OS);
	makeCopySave(histNTracks2OS);
	normaliseHist(histNTracks2OS);

	makeCopySave(histNTracksAll1);
	normaliseHist(histNTracksAll1);
	makeCopySave(histNTracksAll2);
	normaliseHist(histNTracksAll2);

	makeCopySave(histNSoftTracks1);
	normaliseHist(histNSoftTracks1);
	makeCopySave(histNSoftTracks2);
	normaliseHist(histNSoftTracks2);
	makeCopySave(histNSoftTracks1OS);
	normaliseHist(histNSoftTracks1OS);
	makeCopySave(histNSoftTracks2OS);
	normaliseHist(histNSoftTracks2OS);

	// Cumulative plots
	histNTracksCum1   = (TH1D*)histNTracksAbs1->Clone("hNTracksCum1");
	histNTracksCum2   = (TH1D*)histNTracksAbs2->Clone("hNTracksCum2");
	histNTracksCum1OS = (TH1D*)histNTracksAbs1OS->Clone("hNTracksCum1OS");
	histNTracksCum2OS = (TH1D*)histNTracksAbs2OS->Clone("hNTracksCum2OS");

	histNSoftTracksCum1   = (TH1D*)histNSoftTracksAbs1->Clone("hNSoftTracksCum1");
	histNSoftTracksCum2   = (TH1D*)histNSoftTracksAbs2->Clone("hNSoftTracksCum2");
	histNSoftTracksCum1OS = (TH1D*)histNSoftTracksAbs1OS->Clone("hNSoftTracksCum1OS");
	histNSoftTracksCum2OS = (TH1D*)histNSoftTracksAbs2OS->Clone("hNSoftTracksCum2OS");

	for (int i = 1; i <= histNTracks1->GetNbinsX(); i++){
		histNTracksCum1->SetBinContent(i,histNTracksCum1->GetBinContent(i-1) + histNTracksAbs1->GetBinContent(i));
		histNTracksCum2->SetBinContent(i,histNTracksCum2->GetBinContent(i-1) + histNTracksAbs2->GetBinContent(i));
		histNTracksCum1OS->SetBinContent(i,histNTracksCum1OS->GetBinContent(i-1) + histNTracksAbs1OS->GetBinContent(i));
		histNTracksCum2OS->SetBinContent(i,histNTracksCum2OS->GetBinContent(i-1) + histNTracksAbs2OS->GetBinContent(i));

		histNSoftTracksCum1->SetBinContent(i,histNSoftTracksCum1->GetBinContent(i-1) + histNSoftTracksAbs1->GetBinContent(i));
		histNSoftTracksCum2->SetBinContent(i,histNSoftTracksCum2->GetBinContent(i-1) + histNSoftTracksAbs2->GetBinContent(i));
		histNSoftTracksCum1OS->SetBinContent(i,histNSoftTracksCum1OS->GetBinContent(i-1) + histNSoftTracksAbs1OS->GetBinContent(i));
		histNSoftTracksCum2OS->SetBinContent(i,histNSoftTracksCum2OS->GetBinContent(i-1) + histNSoftTracksAbs2OS->GetBinContent(i));
	}

	cout << "# muons with 1+ track with pT >1 GeV: " << n1 << endl;
	cout << "# muons with 1+ track with pT >2.5 GeV: " << n2p5 << endl;
	cout << "# muons with 1+ OS track with pT > 2.5GeV: " << n2p5OS << endl;
	cout << "nMuPass: " << nMuPass << endl;



	// app += "_samePtEta";

	drawHistAndSave(histMu1Pt, "HISTE", "Mu1Pt", directory, app);
	drawHistAndSave(histMu2Pt, "HISTE", "Mu2Pt", directory, app);
	drawHistAndSave(histTrack1Pt, "HISTE", "Track1Pt", directory, app);
	drawHistAndSave(histTrack2Pt, "HISTE", "Track2Pt", directory, app);

	drawHistAndSave(histSys1Pt, "HISTE", "Sys1Pt", directory, app);
	drawHistAndSave(histSys2Pt, "HISTE", "Sys2Pt", directory, app);
	drawHistAndSave(histDRSys, "HISTE", "DRSys", directory, app);
	drawHistAndSave(histDEtaVsDPhiSys, "COLZ", "DEtaVsDPhiSys", directory, app);
	
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

	// All track with pT > 1
	drawHistAndSave(histNTracksAll1, "HISTE", "NTracks1All_NS", directory, app);
	drawHistAndSave(histNTracksAll2, "HISTE", "NTracks2All_NS", directory, app);
	drawHistAndSave(histNTracksAllAbs1, "HISTE", "NTracks1AllAbs_NS", directory, app);
	drawHistAndSave(histNTracksAllAbs2, "HISTE", "NTracks2AllAbs_NS", directory, app);
	
	// Soft tracks around muon
	drawHistAndSave(histNSoftTracks1, "HISTE", "NSoftTracks1_NS", directory, app);
	drawHistAndSave(histNSoftTracks2, "HISTE", "NSoftTracks2_NS", directory, app);
	drawHistAndSave(histNSoftTracks1OS, "HISTE", "NSoftTracks1_OS", directory, app);
	drawHistAndSave(histNSoftTracks2OS, "HISTE", "NSoftTracks2_OS", directory, app);

	drawHistAndSave(histNSoftTracksAbs1, "HISTE", "NSoftTracksAbs1_NS", directory, app);
	drawHistAndSave(histNSoftTracksAbs2, "HISTE", "NSoftTracksAbs2_NS", directory, app);
	drawHistAndSave(histNSoftTracksAbs1OS, "HISTE", "NSoftTracksAbs1_OS", directory, app);
	drawHistAndSave(histNSoftTracksAbs2OS, "HISTE", "NSoftTracksAbs2_OS", directory, app);

	drawHistAndSave(histNSoftTracksCum1, "HISTE", "NSoftTracksCum1_NS", directory, app);
	drawHistAndSave(histNSoftTracksCum2, "HISTE", "NSoftTracksCum2_NS", directory, app);
	drawHistAndSave(histNSoftTracksCum1OS, "HISTE", "NSoftTracksCum1_OS", directory, app);
	drawHistAndSave(histNSoftTracksCum2OS, "HISTE", "NSoftTracksCum2_OS", directory, app);

	drawHistAndSave(histDRMuMu, "HISTE", "DRMuMu", directory, app);
	drawHistAndSave(histDEtaVsDPhiMuMu, "COLZ", "DEtaVsDPhiMuMu", directory, app);

	if (doSignal){
		drawHistAndSave(histSys1PtTruth, "HISTE", "Sys1PtTruth", directory, app);
		drawHistAndSave(histSys2PtTruth, "HISTE", "Sys2PtTruth", directory, app);
		drawHistAndSave(histDRSysTruth, "HISTE", "DRSysTruth", directory, app);
		drawHistAndSave(histDEtaVsDPhiSysTruth, "COLZ", "DEtaVsDPhiSysTruth", directory, app);
		drawHistAndSave(histNuPt, "HISTE", "NuPt", directory, app);
		drawHistAndSave(histDRa1, "HISTE", "DRa1", directory, app);
		drawHistAndSave(histDRa2, "HISTE", "DRa2", directory, app);
		drawHistAndSave(histPID, "HISTE", "PID", directory, app);
		drawHistAndSave(histHPt, "HISTE", "HPt", directory, app);
	}

	TArc problemRing1(0,0,1,0,90);
	problemRing1.SetLineColor(kRed);
	problemRing1.SetLineWidth(2);
	problemRing1.SetFillStyle(0);

	TArc problemRing2(0,0,2,0,90);
	problemRing2.SetLineColor(kRed);
	problemRing2.SetLineWidth(2);
	problemRing2.SetFillStyle(0);

	histTkEtaVsPhi1->Draw("COLZ");
	// problemRing1.Draw("only");
	c.SaveAs((directory+"/TkEtaVsPhi1_"+delph+"_"+app+".pdf").c_str());
	histTkEtaVsPhi2->Draw("COLZ");
	// problemRing1.Draw("only");
	// problemRing2.Draw("only");
	c.SaveAs((directory+"/TkEtaVsPhi2_"+delph+"_"+app+".pdf").c_str());

	// drawHistAndSave(histRand, "HISTE", "RandTest", directory, app);

	///////////////////////////////
	// Write hists to ROOT file //
	///////////////////////////////

	histMu1Pt->Write("",TObject::kOverwrite);
	histMu2Pt->Write("",TObject::kOverwrite);
	histTrack1Pt->Write("",TObject::kOverwrite);
	histTrack2Pt->Write("",TObject::kOverwrite);
	histSys1Pt->Write("",TObject::kOverwrite);
	histSys2Pt->Write("",TObject::kOverwrite);
	histDRSys->Write("",TObject::kOverwrite);
	histDEtaVsDPhiSys->Write("",TObject::kOverwrite);
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
	histNTracksAll1->Write("",TObject::kOverwrite);
	histNTracksAll2->Write("",TObject::kOverwrite);
	histNTracksAllAbs1->Write("",TObject::kOverwrite);
	histNTracksAllAbs2->Write("",TObject::kOverwrite);

	histNSoftTracks1->Write("",TObject::kOverwrite);
	histNSoftTracks2->Write("",TObject::kOverwrite);
	histNSoftTracks1OS->Write("",TObject::kOverwrite);
	histNSoftTracks2OS->Write("",TObject::kOverwrite);
	histNSoftTracksAbs1->Write("",TObject::kOverwrite);
	histNSoftTracksAbs2->Write("",TObject::kOverwrite);
	histNSoftTracksAbs1OS->Write("",TObject::kOverwrite);
	histNSoftTracksAbs2OS->Write("",TObject::kOverwrite);
	histNSoftTracksCum1->Write("",TObject::kOverwrite);
	histNSoftTracksCum2->Write("",TObject::kOverwrite);
	histNSoftTracksCum1OS->Write("",TObject::kOverwrite);
	histNSoftTracksCum2OS->Write("",TObject::kOverwrite);

	histDRMuMu->Write("",TObject::kOverwrite);
	histDEtaVsDPhiMuMu->Write("",TObject::kOverwrite);

	if (doSignal){
		histSys1PtTruth->Write("",TObject::kOverwrite);
		histSys2PtTruth->Write("",TObject::kOverwrite);
		histDRSysTruth->Write("",TObject::kOverwrite);
		histDEtaVsDPhiSysTruth->Write("",TObject::kOverwrite);
		histHPt->Write("",TObject::kOverwrite);
		histNuPt->Write("",TObject::kOverwrite);
		histDRa1->Write("",TObject::kOverwrite);
		histDRa2->Write("",TObject::kOverwrite);
		histPID->Write("",TObject::kOverwrite);
	}

	histTkEtaVsPhi1->Write("",TObject::kOverwrite);
	histTkEtaVsPhi2->Write("",TObject::kOverwrite);
	
	outFile->Close();

}
