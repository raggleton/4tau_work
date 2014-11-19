#include "commonFunctions.h"
#include "classes/SortableObject.h"
#include "cuts.h"
#include "rescale.h"
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

	MCsource source      = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	bool doSignal        = pOpts.getSignal(); // for signal or not
	// bool doMu         = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly  = pOpts.getMuOrdering(); // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT           = pOpts.getHLT(); // whether to use MC that has HLT cuts already applied or not.
	// bool DEBUG        = pOpts.getVerbose(); // output debug statments
	double deltaR        = pOpts.getdR(); // dR(mu-mu) value to use
	bool doRescale       = pOpts.doRescaling(); // whether to rescale signal tracks.
	bool do1to1p5        = false; // for additional sideband studies. Slower?

	// for rescaling signal tracks
	Rescaler r;

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
	// std::vector<double> massBins {0,1,2,3,10};
	// int nBinsX = massBins.size()-1;

	// ------------------------
	// Declare hists
	// ------------------------

	TH1D* histDRmutk_hard                = new TH1D("hDRmutk_hard", "#DeltaR(#mu-tk) (hard)",40,0,2);
	TH1D* histDRmutk_soft                = new TH1D("hDRmutk_soft", "#DeltaR(#mu-tk) (soft)",40,0,2);

	TH1D* histDRmutk_hard_noniso         = new TH1D("hDRmutk_hard_noniso", "#DeltaR(#mu-tk) (hard)",40,0,2);
	TH1D* histDRmutk_soft_noniso         = new TH1D("hDRmutk_soft_noniso", "#DeltaR(#mu-tk) (soft)",40,0,2);

	// NB SHOULD REALLY RENAME THESE HARD OR SOFT, NOT 1 OR 2

	// mu1/2 pt/eta
	TH1F* histMuHardPt_fine_signal       = new TH1F("hMu1Pt_fine_signal","#mu_{hard} p_{T}; #mu_{hard} p_{T} [GeV];A.U.",100,0.,50.);
	TH1F* histMuSoftPt_fine_signal       = new TH1F("hMu2Pt_fine_signal","#mu_{soft} p_{T}; #mu_{soft} p_{T} [GeV];A.U.",100,0.,50.);
	TH1F* histMuHardEta_fine_signal      = new TH1F("hMu1Eta_fine_signal","#mu_{hard} #eta; #mu_{hard} #eta;A.U.",200,-5.,5.);
	TH1F* histMuSoftEta_fine_signal      = new TH1F("hMu2Eta_fine_signal","#mu_{soft} #eta; #mu_{soft} #eta;A.U.",200,-5.,5.);

	// mu1/2 pt/eta
	TH1F* histMuHardPt_fine_IPSSDR       = new TH1F("hMu1Pt_fine_IPSSDR","#mu_{1} p_{T}; #mu_{1} p_{T} [GeV];A.U.",100,0.,50.);
	TH1F* histMuSoftPt_fine_IPSSDR       = new TH1F("hMu2Pt_fine_IPSSDR","#mu_{2} p_{T}; #mu_{2} p_{T} [GeV];A.U.",100,0.,50.);
	TH1F* histMuHardEta_fine_IPSSDR      = new TH1F("hMu1Eta_fine_IPSSDR","#mu_{1} #eta; #mu_{1} #eta;A.U.",200,-5.,5.);
	TH1F* histMuSoftEta_fine_IPSSDR      = new TH1F("hMu2Eta_fine_IPSSDR","#mu_{2} #eta; #mu_{2} #eta;A.U.",200,-5.,5.);

	// mu1 vs mu 2 pT
	TH2F * HardMuonPtSoftMuonPt_DimuonsH = new TH2F("HardMuonPtSoftMuonPt_DimuonsH","",40,0,200,40,0,200);


	// mu1/2 pt/eta
	// TH1D* histTk1Pt_fine_IPSS         = new TH1D("hTk1Pt_fine_IPSS","Tk_{1} p_{T}; Tk_{1} p_{T} [GeV];A.U.",100,0.,100.);
	// TH1D* histTk2Pt_fine_IPSS         = new TH1D("hTk2Pt_fine_IPSS","Tk_{2} p_{T}; Tk_{2} p_{T} [GeV];A.U.",100,0.,100.);
	// TH1D* histTk1Eta_fine_IPSS        = new TH1D("hTk1Eta_fine_IPSS","Tk_{1} #eta; Tk_{1} #eta;A.U.",200,-5.,5.);
	// TH1D* histTk2Eta_fine_IPSS        = new TH1D("hTk2Eta_fine_IPSS","Tk_{2} #eta; Tk_{2} #eta;A.U.",200,-5.,5.);

	// Tk1/2 pt/eta
	TH1D* histTauTk1Pt                   = new TH1D("hTauTk1Pt","Tk_{1} p_{T}; Tk_{1} p_{T} [GeV];A.U.",100,0.,50.);
	TH1D* histTauTk1Eta                  = new TH1D("hTauTk1Eta","Tk_{1} #eta; Tk_{1} #eta;A.U.",200,-5.,5.);
	TH1D* histTauTk2Pt                   = new TH1D("hTauTk2Pt","Tk_{2} p_{T}; Tk_{2} p_{T} [GeV];A.U.",100,0.,50.);
	TH1D* histTauTk2Eta                  = new TH1D("hTauTk2Eta","Tk_{2} #eta; Tk_{2} #eta;A.U.",200,-5.,5.);
	TH2F* histTauTk1vs2Pt                = new TH2F("hTauTk1vs2Pt","",100,0.,50.,100,0.,50.);
	TH1D* histTauTk12Pt                  = new TH1D("hTauTk12Pt","Tk_{1} p_{T}; Tk_{1} p_{T} [GeV];A.U.",100,0.,50.);
	TH1D* histTauTk12Eta                 = new TH1D("hTauTk12Eta","Tk_{2} #eta; Tk_{2} #eta;A.U.",200,-5.,5.);

	// Tk1/2 pt/eta when mu non-iso (2 or 3 tracks)
	TH1D* histNonIsoTk1Pt                = new TH1D("hNonIsoTk1Pt","Tk_{1} p_{T}; Tk_{1} p_{T} [GeV];A.U.",100,0.,50.);
	TH1D* histNonIsoTk1Eta               = new TH1D("hNonIsoTk1Eta","Tk_{1} #eta; Tk_{1} #eta;A.U.",200,-5.,5.);
	TH1D* histNonIsoTk2Pt                = new TH1D("hNonIsoTk2Pt","Tk_{2} p_{T}; Tk_{2} p_{T} [GeV];A.U.",100,0.,50.);
	TH1D* histNonIsoTk2Eta               = new TH1D("hNonIsoTk2Eta","Tk_{2} #eta; Tk_{2} #eta;A.U.",200,-5.,5.);
	TH2F* histNonIsoTk1vs2Pt             = new TH2F("hNonIsoTk1vs2Pt","",100,0.,50.,100,0.,50.);
	TH1D* histNonIsoTk12Pt               = new TH1D("hNonIsoTk12Pt","Tk_{1} p_{T}; Tk_{1} p_{T} [GeV];A.U.",100,0.,50.);
	TH1D* histNonIsoTk12Eta              = new TH1D("hNonIsoTk12Eta","Tk_{2} #eta; Tk_{2} #eta;A.U.",200,-5.,5.);



	// int nMu(0);
	// int n1(0), n2(0), nMuPass(0);

	///////////////////////
	// Loop over events  //
	///////////////////////
	Long64_t numberOfEntries = getNumberEvents(treeReader, &pOpts);
	cout << "Running over " << numberOfEntries << " events" << endl;

	bool stop = false; // used to stop the loop, for debugging/testing
	int last_pc = -1;
	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// output something every 5%
		int pc = int(100*entry/double(numberOfEntries));
		if(pc%5 == 0 && pc != last_pc) {
			std::cout << pc << "% progress (event " << entry << ")" << endl;
			last_pc = pc;
		}

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
					// 	histMuHardPt_truth_0to1->Fill(muTruth1->PT);
					// } else if (m2 < 2.){
					// 	histM1_truth_1to2->Fill(m1);
					// 	histMuHardPt_truth_1to2->Fill(muTruth1->PT);
					// } else if (m2 < 3.){
					// 	histM1_truth_2to3->Fill(m1);
					// 	histMuHardPt_truth_2to3->Fill(muTruth1->PT);
					// } else{
					// 	histM1_truth_3toInf->Fill(m1);
					// 	histMuHardPt_truth_3toInf->Fill(muTruth1->PT);
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


		// For rest of plots we want our muons to pass signal selection
		// Make pairs, see if they pass all cuts (SS, eta, deltaR, dZ, d0)
		// If they do, store in mu1 and mu2 (mu1 has higher pT)
		std::pair<Track*, Track*> p = testMuons(muons17toInf, muons10to17, &checkMuons, deltaR);
		Track* mu1 = p.first;
		Track* mu2 = p.second;

		if (!(p.first && p.second)) continue;

		histMuHardPt_fine_IPSSDR->Fill(p.first->PT);
		histMuHardEta_fine_IPSSDR->Fill(p.first->Eta);
		histMuSoftPt_fine_IPSSDR->Fill(p.second->PT);
		histMuSoftEta_fine_IPSSDR->Fill(p.second->Eta);
		HardMuonPtSoftMuonPt_DimuonsH->Fill(p.first->PT, p.second->PT);

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


				if(doRescale) {
					r.rescaleTrack(tk1_2p5_OS.at(0), mu1);
					r.rescaleTrack(tk2_2p5_OS.at(0), mu2);
				}

				TLorentzVector track1Mom=tk1_2p5_OS[0]->P4();
				TLorentzVector track2Mom=tk2_2p5_OS[0]->P4();

				// random since mu1Mom andmu2Mom are randomly assigned (if selected at top)
				// double m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
				// double m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

				if (swapped){
					histMuHardPt_fine_signal->Fill(mu2->PT);
					histMuHardEta_fine_signal->Fill(mu2->Eta);
					histMuSoftPt_fine_signal->Fill(mu1->PT);
					histMuSoftEta_fine_signal->Fill(mu1->Eta);

					histTauTk1Pt->Fill(tk2_2p5_OS[0]->PT);
					histTauTk1Eta->Fill(tk2_2p5_OS[0]->Eta);
					histTauTk2Pt->Fill(tk1_2p5_OS[0]->PT);
					histTauTk2Eta->Fill(tk1_2p5_OS[0]->Eta);
					// histTauTk1vs2Pt->Fill(tk2_2p5_OS[0]->PT, tk1_2p5_OS[0]->PT);
					histTauTk12Pt->Fill(tk1_2p5_OS[0]->PT),
					histTauTk12Pt->Fill(tk2_2p5_OS[0]->PT);
					histTauTk12Eta->Fill(tk1_2p5_OS[0]->Eta);
					histTauTk12Eta->Fill(tk2_2p5_OS[0]->Eta);
					histDRmutk_hard->Fill(track2Mom.DeltaR(mu2Mom));
					histDRmutk_soft->Fill(track1Mom.DeltaR(mu1Mom));
				} else {
					histMuHardPt_fine_signal->Fill(mu1->PT);
					histMuHardEta_fine_signal->Fill(mu1->Eta);
					histMuSoftPt_fine_signal->Fill(mu2->PT);
					histMuSoftEta_fine_signal->Fill(mu2->Eta);

					histTauTk1Pt->Fill(tk1_2p5_OS[0]->PT);
					histTauTk1Eta->Fill(tk1_2p5_OS[0]->Eta);
					histTauTk2Pt->Fill(tk2_2p5_OS[0]->PT);
					histTauTk2Eta->Fill(tk2_2p5_OS[0]->Eta);
					// histTauTk1vs2Pt->Fill(tk1_2p5_OS[0]->PT, tk2_2p5_OS[0]->PT);
					histTauTk12Pt->Fill(tk1_2p5_OS[0]->PT),
					histTauTk12Pt->Fill(tk2_2p5_OS[0]->PT);
					histTauTk12Eta->Fill(tk1_2p5_OS[0]->Eta);
					histTauTk12Eta->Fill(tk2_2p5_OS[0]->Eta);
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

			// double m1(0);
			TLorentzVector track1Mom=tk1_2p5_OS[0]->P4();
			TLorentzVector track2Mom=tk2_2p5_OS[0]->P4();

			if ( //tk2_1.size() == 2
				tk2_1to2p5_taucand.size() == 1
				&& tk2_1to2p5.size() == 1
				) {
					// m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();


				if (swapped){
					histNonIsoTk1Pt->Fill(tk2_2p5_OS[0]->PT);
					histNonIsoTk1Eta->Fill(tk2_2p5_OS[0]->Eta);
					histNonIsoTk2Pt->Fill(tk1_2p5_OS[0]->PT);
					histNonIsoTk2Eta->Fill(tk1_2p5_OS[0]->Eta);
					// histNonIsoTk1vs2Pt->Fill(tk2_2p5_OS[0]->PT, tk1_2p5_OS[0]->PT);
					histNonIsoTk12Pt->Fill(tk1_2p5_OS[0]->PT),
					histNonIsoTk12Pt->Fill(tk2_2p5_OS[0]->PT);
					histNonIsoTk12Eta->Fill(tk1_2p5_OS[0]->Eta);
					histNonIsoTk12Eta->Fill(tk2_2p5_OS[0]->Eta);
					histDRmutk_hard_noniso->Fill(track2Mom.DeltaR(mu2Mom));
					histDRmutk_soft_noniso->Fill(track1Mom.DeltaR(mu1Mom));
				} else {
					histNonIsoTk1Pt->Fill(tk1_2p5_OS[0]->PT);
					histNonIsoTk1Eta->Fill(tk1_2p5_OS[0]->Eta);
					histNonIsoTk2Pt->Fill(tk2_2p5_OS[0]->PT);
					histNonIsoTk2Eta->Fill(tk2_2p5_OS[0]->Eta);
					// histNonIsoTk1vs2Pt->Fill(tk1_2p5_OS[0]->PT, tk2_2p5_OS[0]->PT);
					histNonIsoTk12Pt->Fill(tk1_2p5_OS[0]->PT),
					histNonIsoTk12Pt->Fill(tk2_2p5_OS[0]->PT);
					histNonIsoTk12Eta->Fill(tk1_2p5_OS[0]->Eta);
					histNonIsoTk12Eta->Fill(tk2_2p5_OS[0]->Eta);
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
						histNonIsoTk1Pt->Fill(tk2_2p5_OS[0]->PT);
						histNonIsoTk1Eta->Fill(tk2_2p5_OS[0]->Eta);
						histNonIsoTk2Pt->Fill(tk1_2p5_OS[0]->PT);
						histNonIsoTk2Eta->Fill(tk1_2p5_OS[0]->Eta);
						// histNonIsoTk1vs2Pt->Fill(tk2_2p5_OS[0]->PT, tk1_2p5_OS[0]->PT);
						histNonIsoTk12Pt->Fill(tk1_2p5_OS[0]->PT),
						histNonIsoTk12Pt->Fill(tk2_2p5_OS[0]->PT);
						histNonIsoTk12Eta->Fill(tk1_2p5_OS[0]->Eta);
						histNonIsoTk12Eta->Fill(tk2_2p5_OS[0]->Eta);
						histDRmutk_hard_noniso->Fill(track2Mom.DeltaR(mu2Mom));
						histDRmutk_soft_noniso->Fill(track1Mom.DeltaR(mu1Mom));
					} else {
						histNonIsoTk1Pt->Fill(tk1_2p5_OS[0]->PT);
						histNonIsoTk1Eta->Fill(tk1_2p5_OS[0]->Eta);
						histNonIsoTk2Pt->Fill(tk2_2p5_OS[0]->PT);
						histNonIsoTk2Eta->Fill(tk2_2p5_OS[0]->Eta);
						// histNonIsoTk1vs2Pt->Fill(tk1_2p5_OS[0]->PT, tk2_2p5_OS[0]->PT);
						histNonIsoTk12Pt->Fill(tk1_2p5_OS[0]->PT),
						histNonIsoTk12Pt->Fill(tk2_2p5_OS[0]->PT);
						histNonIsoTk12Eta->Fill(tk1_2p5_OS[0]->Eta);
						histNonIsoTk12Eta->Fill(tk2_2p5_OS[0]->Eta);
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


	//////////////////////////////////////////////
	// Print out integrals before we normalise //
	//////////////////////////////////////////////
	// printIntegral(histMuHardPt_fine_IPSS);
	// printIntegral(histMuSoftPt_fine_IPSS);
	// printIntegral(histMuHardEta_fine_IPSS);
	// printIntegral(histMuSoftEta_fine_IPSS);

	printIntegral(histMuHardPt_fine_IPSSDR);
	printIntegral(histMuSoftPt_fine_IPSSDR);
	printIntegral(histMuHardEta_fine_IPSSDR);
	printIntegral(histMuSoftEta_fine_IPSSDR);

	printIntegral(histMuHardPt_fine_signal);
	printIntegral(histMuSoftPt_fine_signal);
	printIntegral(histMuHardEta_fine_signal);
	printIntegral(histMuSoftEta_fine_signal);

	printIntegral(histTauTk1Pt);
	printIntegral(histTauTk1Eta);
	printIntegral(histTauTk2Pt);
	printIntegral(histTauTk2Eta);
	printIntegral(histTauTk12Pt);
	printIntegral(histTauTk12Eta);
	printIntegral(histNonIsoTk1Pt);
	printIntegral(histNonIsoTk1Eta);
	printIntegral(histNonIsoTk2Pt);
	printIntegral(histNonIsoTk2Eta);
	printIntegral(histNonIsoTk12Pt);
	printIntegral(histNonIsoTk12Eta);
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

	if (doRescale) {
		app += "_rescaleQuantile";
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

	/////////////////////////////////
	// PLOT AND SAVE HISTS TO PDF //
	/////////////////////////////////

	// mu1,2 pt,eta
	// drawHistAndSave(histMuHardPt_fine_IPSS, "HISTE", "histMuHardPt_fine_IPSS", directory, app);
	// drawHistAndSave(histMuHardEta_fine_IPSS, "HISTE", "histMuHardEta_fine_IPSS", directory, app);
	// drawHistAndSave(histMuSoftPt_fine_IPSS, "HISTE", "histMuSoftPt_fine_IPSS", directory, app);
	// drawHistAndSave(histMuSoftEta_fine_IPSS, "HISTE", "histMuSoftEta_fine_IPSS", directory, app);
	drawHistAndSave(histMuHardPt_fine_IPSSDR, "HISTE", "histMuHardPt_fine_IPSSDR", directory, app);
	drawHistAndSave(histMuHardEta_fine_IPSSDR, "HISTE", "histMuHardEta_fine_IPSSDR", directory, app);
	drawHistAndSave(histMuSoftPt_fine_IPSSDR, "HISTE", "histMuSoftPt_fine_IPSSDR", directory, app);
	drawHistAndSave(histMuSoftEta_fine_IPSSDR, "HISTE", "histMuSoftEta_fine_IPSSDR", directory, app);

	drawHistAndSave(histMuHardPt_fine_signal, "HISTE", "histMuHardPt_fine_signal", directory, app);
	drawHistAndSave(histMuHardEta_fine_signal, "HISTE", "histMuHardEta_fine_signal", directory, app);
	drawHistAndSave(histMuSoftPt_fine_signal, "HISTE", "histMuSoftPt_fine_signal", directory, app);
	drawHistAndSave(histMuSoftEta_fine_signal, "HISTE", "histMuSoftEta_fine_signal", directory, app);

	drawHistAndSave(HardMuonPtSoftMuonPt_DimuonsH, "COLZ", "HardMuonPtSoftMuonPt_DimuonsH", directory, app);

	// tk1,2 pt, eta
	drawHistAndSave(histTauTk1Pt, "HISTE", "histTauTk1Pt", directory, app);
	drawHistAndSave(histTauTk1Eta, "HISTE", "histTauTk1Eta", directory, app);
	drawHistAndSave(histTauTk2Pt, "HISTE", "histTauTk2Pt", directory, app);
	drawHistAndSave(histTauTk2Eta, "HISTE", "histTauTk2Eta", directory, app);
	drawHistAndSave(histTauTk1vs2Pt, "COLZ", "histTauTk1vs2Pt", directory, app);
	drawHistAndSave(histTauTk12Pt, "HISTE", "histTauTk12Pt", directory, app);
	drawHistAndSave(histTauTk12Eta, "HISTE", "histTauTk12Eta", directory, app);

	drawHistAndSave(histNonIsoTk1Pt, "HISTE", "histNonIsoTk1Pt", directory, app);
	drawHistAndSave(histNonIsoTk1Eta, "HISTE", "histNonIsoTk1Eta", directory, app);
	drawHistAndSave(histNonIsoTk2Pt, "HISTE", "histNonIsoTk2Pt", directory, app);
	drawHistAndSave(histNonIsoTk2Eta, "HISTE", "histNonIsoTk2Eta", directory, app);
	drawHistAndSave(histNonIsoTk1vs2Pt, "COLZ", "histNonIsoTk12Pt", directory, app);
	drawHistAndSave(histNonIsoTk12Pt, "HISTE", "histNonIsoTk12Pt", directory, app);
	drawHistAndSave(histNonIsoTk12Eta, "HISTE", "histNonIsoTk12Eta", directory, app);

	drawHistAndSave(histDRmutk_hard, "HISTE", "DRmutk_hard", directory, app);
	drawHistAndSave(histDRmutk_soft, "HISTE", "DRmutk_soft", directory, app);
	drawHistAndSave(histDRmutk_hard_noniso, "HISTE", "DRmutk_hard_noniso", directory, app);
	drawHistAndSave(histDRmutk_soft_noniso, "HISTE", "DRmutk_soft_noniso", directory, app);

	//////////////////////////
	// Write hists to file //
	//////////////////////////

	histMuHardPt_fine_IPSSDR->Write("", TObject::kOverwrite);
	histMuHardEta_fine_IPSSDR->Write("", TObject::kOverwrite);
	histMuSoftPt_fine_IPSSDR->Write("", TObject::kOverwrite);
	histMuSoftEta_fine_IPSSDR->Write("", TObject::kOverwrite);

	histMuHardPt_fine_signal->Write("", TObject::kOverwrite);
	histMuHardEta_fine_signal->Write("", TObject::kOverwrite);
	histMuSoftPt_fine_signal->Write("", TObject::kOverwrite);
	histMuSoftEta_fine_signal->Write("", TObject::kOverwrite);

	HardMuonPtSoftMuonPt_DimuonsH->Write("", TObject::kOverwrite);

	histTauTk1Pt->Write("", TObject::kOverwrite);
	histTauTk1Eta->Write("", TObject::kOverwrite);
	histTauTk2Pt->Write("", TObject::kOverwrite);
	histTauTk2Eta->Write("", TObject::kOverwrite);
	histTauTk1vs2Pt->Write("", TObject::kOverwrite);
	histTauTk12Pt->Write("", TObject::kOverwrite);
	histTauTk12Eta->Write("", TObject::kOverwrite);

	histNonIsoTk1Pt->Write("", TObject::kOverwrite);
	histNonIsoTk1Eta->Write("", TObject::kOverwrite);
	histNonIsoTk2Pt->Write("", TObject::kOverwrite);
	histNonIsoTk2Eta->Write("", TObject::kOverwrite);
	histNonIsoTk1vs2Pt->Write("", TObject::kOverwrite);
	histNonIsoTk12Pt->Write("", TObject::kOverwrite);
	histNonIsoTk12Eta->Write("", TObject::kOverwrite);

	histDRmutk_hard->Write("",TObject::kOverwrite);
	histDRmutk_soft->Write("",TObject::kOverwrite);
	histDRmutk_hard_noniso->Write("",TObject::kOverwrite);
	histDRmutk_soft_noniso->Write("",TObject::kOverwrite);

	if (doSignal){
	}

	outFile->Close();

	delete treeReader;
}
