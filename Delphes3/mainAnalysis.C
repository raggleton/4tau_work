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
void mainAnalysis(int argc, char* argv[])
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
	
	// Track distributions
	TH1D *histNTracks1           = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks2           = new TH1D("hNTracks2" ,"Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);
	TH1D *histNTracks1OS         = new TH1D("hNTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	TH1D *histNTracks2OS         = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);

	TH1D *histNTracksCum1        = new TH1D("hNTracksCum1" ,"Cumulative Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNTracksCum2        = new TH1D("hNTracksCum2" ,"Cumulative Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);
	TH1D *histNTracksCum1OS      = new TH1D("hNTracksCum1OS" ,"Cumulative Number of tracks about mu1, OS,p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); Average N_{trk} about #mu_{1}", 50,0,5);
	TH1D *histNTracksCum2OS      = new TH1D("hNTracksCum2OS" ,"Cumulative Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); Average N_{trk} about #mu_{2}", 50,0,5);

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

	// Plots for testing invariant mass correlation
	double massBins[6]           = {0,1,2,3,4,10};
	TH1D* histM1                 = new TH1D("hM1", "Inv. Mass of 1st system, full selection; m(#mu_{1}-tk) [GeV]; N_{events}",10,0,10);
	TH1D* histM2                 = new TH1D("hM2", "Inv. Mass of 2st system, full selection; m(#mu_{2}-tk) [GeV]; N_{events}",10,0,10);

	// MC truth - use actual mu-tk pairs from tau
	TH1D* histM1_truth_0to1      = new TH1D("hM1_truth_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_truth_1to2      = new TH1D("hM1_truth_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_truth_2to3      = new TH1D("hM1_truth_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_truth_3toInf    = new TH1D("hM1_truth_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	
	// actual dist using signal selection
	TH1D* histM1_0to1            = new TH1D("hM1_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_1to2            = new TH1D("hM1_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_2to3            = new TH1D("hM1_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_3toInf          = new TH1D("hM1_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);

	// actual dist using sideband selection
	TH1D* histM1_side_0to1       = new TH1D("hM1_side_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_side_1to2       = new TH1D("hM1_side_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_side_2to3       = new TH1D("hM1_side_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);
	TH1D* histM1_side_3toInf     = new TH1D("hM1_side_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",5,massBins);

	// mu1 pT plots in bins of m2
	TH1D* histMu1Pt_0to1         = new TH1D("hMu1Pt_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_1to2         = new TH1D("hMu1Pt_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_2to3         = new TH1D("hMu1Pt_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_3toInf       = new TH1D("hMu1Pt_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);

	// mu1 pT plots in bins of m2 - MC truth
	TH1D* histMu1Pt_truth_0to1   = new TH1D("hMu1Pt_truth_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_1to2   = new TH1D("hMu1Pt_truth_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_2to3   = new TH1D("hMu1Pt_truth_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_3toInf = new TH1D("hMu1Pt_truth_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);

	// mu1 pT plots in bins of m2 - sideband
	TH1D* histMu1Pt_side_0to1    = new TH1D("hMu1Pt_side_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_side_1to2    = new TH1D("hMu1Pt_side_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_side_2to3    = new TH1D("hMu1Pt_side_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_side_3toInf  = new TH1D("hMu1Pt_side_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);

	// 2D plots of m1 Vs m2 - sideband
	TH2D* histM1vsM2_side        = new TH2D("hM1vsM2_side","m(#mu_{1}-tk) vs m(#mu_{2}-tk);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",5,massBins,5,massBins);
	TH2D* histM1timesM1_side     = new TH2D("hM1timesM2_side","m(sideband) #times m(sideband);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",5,massBins,5,massBins);

	int nMu(0);
	int n2p5(0), n2p5OS(0); // count # muons with 1+ tracks with pT > 2.5 for SS+OS, and OS
	int nOnly2p5(0), nOnly2p5OS(0); // count # muons with 1 tracks with  pT > 2.5 for SS+OS, and OS
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

		//////////////////////////////////////////////////////////////////////
		// First, get the two highest pT muons in the event, store their pT //
		// and pointers to the GenParticles                                 //
		//////////////////////////////////////////////////////////////////////
		
		GenParticle *cand(nullptr),*mu1(nullptr), *mu2(nullptr);
		// Track *candTk(nullptr);

		// Store pT of highgest and 2nd highest pT muons
		double muLeadingPT(0.);
		double muSubLeadingPT(0.);
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
			// histRand->Fill(randNum);
			if (randNum > 0.5){
				mu1 = origMu2;
				mu2 = origMu1;
			}
		}

		TLorentzVector mu1Mom, mu2Mom;
		mu1Mom = mu1->P4();
		mu2Mom = mu2->P4();

		histMu1Pt->Fill(mu1->PT);
		histMu2Pt->Fill(mu2->PT);

		//////////////////////////////////////////////////////
		// Get the hard interaction particles for signal MC //
		// No selection cuts applied (only >=2 muons)       //
		//////////////////////////////////////////////////////
		
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
			bool atLeastTk2p5 = false; // to monitor if theres a tk with pT > 2.5
			bool atLeastTk2p5OS = false; // same but for OS tk-muon
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

						histNTracks1->Fill(dR1);
						histNTracks2->Fill(dR2);
						atLeastTk2p5 = true;

						histTkEtaVsPhi1->Fill(fabs(candTk->Eta - mu1Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu1Mom)));
						histTkEtaVsPhi2->Fill(fabs(candTk->Eta - mu2Mom.Eta()),fabs((candTk->P4()).DeltaPhi(mu2Mom)));

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
						tk1_1to2p5_alldR.push_back(candTk);
						tk2_1to2p5_alldR.push_back(candTk);
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
			// ANOTHER SIDEBAND REGION
			// where at least one muon has an additional track with 1< pT < 2.5,
			// within dR < 0.5. No sign requirement.
			if (tk1_1.size() >= 1 && tk2_1.size() >= 1 // at least 1 tk with pT > 1 about each muon
				&& tk1_2p5.size() == 1 && tk2_2p5.size() == 1 // main tk with pT > 2.5
				&& (   (tk1_1to2p5.size() == 1 && tk2_1to2p5.size() == 1) // additional soft track about mu1 and/or mu2
					|| (tk1_1to2p5.size() == 0 && tk2_1to2p5.size() == 1) 
					|| (tk1_1to2p5.size() == 1 && tk2_1to2p5.size() == 0))
				){ 
				
				double m1(0), m2(0);				
				if (tk1_1to2p5.size() == 1)
					m1 = (mu1Mom+tk1_2p5[0]->P4()+tk1_1to2p5[0]->P4()).M();
				else
					m1 = (mu1Mom+tk1_2p5[0]->P4()).M();
				
				if(tk2_1to2p5.size() == 1)
					m2 = (mu2Mom+tk2_2p5[0]->P4()+tk2_1to2p5[0]->P4()).M();
				else
					m2 = (mu2Mom+tk2_2p5[0]->P4()).M();
				
				// Fill 2D m1 vs m2 plot
				histM1vsM2_side->Fill(m1,m2);

				// Fill for the various m2 bins
				if(m2 < 1.) {
					histM1_side_0to1->Fill(m1);
					histMu1Pt_side_0to1->Fill(mu1->PT);
				} else if (m2 < 2.) {
					histM1_side_1to2->Fill(m1);
					histMu1Pt_side_1to2->Fill(mu1->PT);
				} else if (m2 < 3.) {
					histM1_side_2to3->Fill(m1);
					histMu1Pt_side_2to3->Fill(mu1->PT);
				} else {
					histM1_side_3toInf->Fill(m1);
					histMu1Pt_side_3toInf->Fill(mu1->PT);
				}

			}

			// Slightly different region - for additional track investigations
			// For soft track distributions
			// ATM it uses signal region. Not put in the signal region bit above, as subject to future modification
			if (tk1_1.size() == 1 && tk2_1.size() == 1 && tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1){
				
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

		} // end of muon selection
		
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

	histNSoftTracksAbs1->Scale(1./nOnly2p5OS);
	histNSoftTracksAbs1OS->Scale(1./nOnly2p5OS);
	histNSoftTracksAbs2->Scale(1./nOnly2p5OS);
	histNSoftTracksAbs2OS->Scale(1./nOnly2p5OS);

	histNSoftTracksAbs1->SetYTitle("Ave. N_{trk} per #mu_{1}");
	histNSoftTracksAbs2->SetYTitle("Ave. N_{trk} per #mu_{2}");
	histNSoftTracksAbs1OS->SetYTitle("Ave. OS N_{trk} per #mu_{1}");
	histNSoftTracksAbs2OS->SetYTitle("Ave. OS N_{trk} per #mu_{2}");

	// AU scaling
	normaliseHist(histNTracks1);
	normaliseHist(histNTracks2);
	normaliseHist(histNTracks1OS);
	normaliseHist(histNTracks2OS);

	normaliseHist(histNSoftTracks1);
	normaliseHist(histNSoftTracks2);
	normaliseHist(histNSoftTracks1OS);
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

	cout << "# muons with 1+ track with pT >2.5 GeV: " << n2p5 << endl;
	cout << "# muons with 1+ OS track with pT > 2.5GeV: " << n2p5OS << endl;
	cout << "nMuPass: " << nMuPass << endl;

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

	// Get directory that input file was in - put plots in there
	std::string directory = getDirectory(chain.GetFile());
	// Get Delphes file config used - last part of directory name
	std::string delph = getDelph(directory);


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

	// Draw mass plots
	// Make m1 sideband plot
	// TH1D* histM1_side = new TH1D("hM1_side","m(#mu_{1}-tk) in sideband;m(#mu_{1}-tk) [GeV];A.U.",5,massBins);
	// histM1_side->Add(histM1_side_0to1);
	// histM1_side->Add(histM1_side_1to2);
	// histM1_side->Add(histM1_side_2to3);
	// histM1_side->Add(histM1_side_3toInf);

	// Do some normalizing
	// normaliseHist(histM1_side);
	// normaliseHist(histM1vsM2_side);
	
	// Do corrections plot by making m1*m1 first, then dividing m1vsm2 by m1*m1
	// Don't need to normalise m1timesm1, as histM1_side already normalised
	// for(int a = 1; a <= 5; a++){
		// for (int b = 1; b <=5; b++){
			// histM1timesM1_side->SetBinContent(a,b,histM1_side->GetBinContent(a)*histM1_side->GetBinContent(b));
		// }
	// }
	// TH2D* histM1vsM2_corrections_side = (TH2D*)histM1vsM2_side->Clone("hM1vsM2_corrections_side");
	// histM1vsM2_corrections_side->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	// histM1vsM2_corrections_side->Divide(histM1timesM1_side);

	// drawHistAndSave(histM1, "HISTE", "M1", directory, app);
	// drawHistAndSave(histM2, "HISTE", "M2", directory, app);

	// drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - signal region;m(#mu_{1}-tk) [GeV]; A.U.", histM1_0to1, histM1_1to2, histM1_2to3, histM1_3toInf, "M1_M2", directory, app);
	// drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - sideband region (at least 1 #mu has add. tk with p_{T} = (1,2.5));m(#mu_{1}-tk) [GeV]; A.U.", histM1_side_0to1, histM1_side_1to2, histM1_side_2to3, histM1_side_3toInf, "M1_M2_side", directory, app);
	if(doSignal){
		// drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - MC truth;m(#mu_{1}-tk) [GeV]; A.U.", histM1_truth_0to1, histM1_truth_1to2, histM1_truth_2to3, histM1_truth_3toInf, "M1_M2_truth", directory, app);
	}

	// drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - signal region;#mu_{1} p_{T} [GeV]; A.U.", histMu1Pt_0to1, histMu1Pt_1to2, histMu1Pt_2to3, histMu1Pt_3toInf, "Mu1Pt_M2", directory, app);
	// drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - MC truth;#mu_{1} p_{T} [GeV]; A.U.", histMu1Pt_truth_0to1, histMu1Pt_truth_1to2, histMu1Pt_truth_2to3, histMu1Pt_truth_3toInf, "Mu1Pt_M2_truth", directory, app);
	// drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - sideband;#mu_{1} p_{T} [GeV]; A.U.", histMu1Pt_side_0to1, histMu1Pt_side_1to2, histMu1Pt_side_2to3, histMu1Pt_side_3toInf, "Mu1Pt_M2_side", directory, app);

	// drawHistAndSave(histM1_side, "HISTE", "M1_side", directory, app);
	// drawHistAndSave(histM1vsM2_side, "colz","M1vsM2_side", directory, app);
	// drawHistAndSave(histM1timesM1_side, "colz","M1timesM1_side", directory, app);
	// drawHistAndSave(histM1vsM2_corrections_side, "colzTEXTE","M1vsM2_corrections_side", directory, app);


	///////////////////////////////
	// Write hists to ROOT file //
	///////////////////////////////
	TFile* outFile = TFile::Open((directory+"/output_"+delph+"_"+app+".root").c_str(),"UPDATE");

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
	
	// Mass plots
	// histM1->Write("",TObject::kOverwrite);
	// histM2->Write("",TObject::kOverwrite);
	// if (doSignal){
	// 	histM1_truth_0to1->Write("",TObject::kOverwrite);
	// 	histM1_truth_1to2->Write("",TObject::kOverwrite);
	// 	histM1_truth_2to3->Write("",TObject::kOverwrite);
	// 	histM1_truth_3toInf->Write("",TObject::kOverwrite);
	// }
	// histM1_0to1->Write("",TObject::kOverwrite);
	// histM1_1to2->Write("",TObject::kOverwrite);
	// histM1_2to3->Write("",TObject::kOverwrite);
	// histM1_3toInf->Write("",TObject::kOverwrite);
	// histM1_side_0to1->Write("",TObject::kOverwrite);
	// histM1_side_1to2->Write("",TObject::kOverwrite);
	// histM1_side_2to3->Write("",TObject::kOverwrite);
	// histM1_side_3toInf->Write("",TObject::kOverwrite);

	outFile->Close();

}
