#include "commonFunctions.h"
#include <boost/program_options.hpp>

using std::cout;
using std::endl;

namespace po = boost::program_options;

/**
 * Template script for Delphes analysis. Use with makeScript.sh
 */
void basicScript(int argc, char* argv[])
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	bool doSignal = true;
	bool doMu = true; // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = true; // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT = true; // whether to use MC that has HLT cuts already applied or not.
	
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("doSignal", po::value<bool>(&doSignal), "TRUE - do signal, FALSE - do QCDb_mu")
		("swapMuRandomly", po::value<bool>(&swapMuRandomly), "TRUE - mu 1,2 randomly assigned, FALSE - mu 1,2 pT ordered")
		("doHLT", po::value<bool>(&doHLT), "TRUE - use samples with HLT_Mu17_Mu8 during generation, FALSE - no HLT cuts")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
	    cout << desc << "\n";
	    return;
	}

	if (vm.count("doSignal")) {
	    doSignal = vm["doSignal"].as<bool>();
	} else {
	    cout << "Signal/QCD was not set. Defaulting to signal.\n";
	}
	if (vm.count("swapMuRandomly")) {
	    swapMuRandomly = vm["swapMuRandomly"].as<bool>();
	} else {
	    cout << "Mu ordering not set. Defaulting to random.\n";
	}
	if (vm.count("doHLT")) {
	    doHLT = vm["doHLT"].as<bool>();
	} else {
	    cout << "HLT requirement not set. Defaulting to using samples with HLT cuts.\n";
	}
	
	
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
					// cout << m1 << "     " << m2 << endl;
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

		double muLeadingPT = 0.;
		double muSubLeadingPT = 0.;
		// Get highest pT muon
		for (int i = 0; i < branchGenMuons->GetEntries(); i++){
			cand = (GenParticle*) branchGenMuons->At(i);
			if (cand->PT > muLeadingPT) {
				mu1 = cand;
				muLeadingPT = cand->PT;
			}
		}
		// Get 2nd highest pT muon
		for(int j = 0; j < branchGenMuons->GetEntries(); j++){
			cand = (GenParticle*) branchGenMuons->At(j);
			if ((cand->PT > muSubLeadingPT) && (cand->PT != mu1->PT)) {
				mu2 = cand;
				muSubLeadingPT = cand->PT;
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
		
		if ((muLeadingPT < 17.)
		|| (muSubLeadingPT < 10.)
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
		// so tk1 is the track nearest to muon1, (may or may not be highest pT, depends if random swapping is on)
		std::vector<Track*> tk1_1;
		std::vector<Track*> tk2_1;

		// same but with pT >2.5
		std::vector<Track*> tk1_2p5;
		std::vector<Track*> tk2_2p5;

		// same but with pT > 2.5, OS to muon
		std::vector<Track*> tk1_2p5_OS;
		std::vector<Track*> tk2_2p5_OS;

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
				}
			} // End of track selection
		} // End of track loop

		// Now pT order the track collections
		// std::sort(tk1.begin(), tk1.end(), sortTracksByPT);
		// std::sort(tk2.begin(), tk2.end(), sortTracksByPT);

		// Can now deal with signal or sideband regions
		
		// SIGNAL SELECTION
		if (tk1_1.size() == 1 && tk2_1.size() == 1 
		&& tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {
			// do something in signal region...
		}		

		// Clean up
		tk1_1.clear();
		tk2_1.clear();
		tk1_2p5.clear();
		tk2_2p5.clear();
		tk1_2p5_OS.clear();
		tk2_2p5_OS.clear();
	} // end of event loop


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

	// app += "_samePtEta";

	// Get directory that input file was in - put plots in there
	std::string directory = getDirectory(chain.GetFile());
	// Get Delphes file config used - last part of directory name
	std::string delph = getDelph(directory);
	
	// drawHistAndSave(histmu1PT, "HISTE", "muLeadingPT", directory, app);
	// drawHistAndSave(histmu2PT, "HISTE", "muSubLeadingPT", directory, app);
	// drawHistAndSave(histTrack1Pt, "HISTE", "Track1Pt", directory, app);
	// drawHistAndSave(histTrack2Pt, "HISTE", "Track2Pt", directory, app);

	// TFile* outFile = TFile::Open((directory+"/output_"+delph+"_"+app+".root").c_str(),"UPDATE");

	// histNMu->Write("",TObject::kOverwrite);
	// histmu1PT->Write("",TObject::kOverwrite);
	// histmu2PT->Write("",TObject::kOverwrite);
	// histmuLeadingPTSel
	// if (doSignal){
	// 	histDRa1->Write("",TObject::kOverwrite);
	// 	histDRa2->Write("",TObject::kOverwrite);
	// 	histPID->Write("",TObject::kOverwrite);
	// }

	// outFile->Close();

}
