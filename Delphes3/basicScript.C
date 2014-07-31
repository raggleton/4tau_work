#include "commonFunctions.h"
#include "cuts.h"
// #include <boost/program_options.hpp>

using std::cout;
using std::endl;


/**
 * Template script for Delphes analysis. Use with makeScript.sh
 */
void basicScript(int argc, char* argv[])
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	ProgramOpts pOpts(argc, argv);

	MCsource source     = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	bool doSignal       = pOpts.getSignal(); // true if doing signal
	bool doMu           = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
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
	// TClonesArray *branchGenMuons = treeReader->UseBranch("OnlyGenMuons");
	TClonesArray *branchGenMuons = treeReader->UseBranch("GenMuon");
	// TClonesArray *branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray *branchAll      = treeReader->UseBranch("AllParticle");

	// Book histograms here
	// TH1D *histNTracks1         = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); A.U.", 50,0,5);
	// TH1D *histNTracks2OS       = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); A.U.", 50,0,5);

	///////////////////////
	// Loop over events //
	///////////////////////
	Long64_t numberOfEntries = getNumberEvents(treeReader, &pOpts);
	cout << "Running over " << numberOfEntries << " events" << endl;
	
	bool stop = false;// used to stop the loop, for debugging/testing

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
			GenParticle *a1(nullptr), *a2(nullptr);
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
			GenParticle *tau1a(nullptr), *tau1b(nullptr), *tau2a(nullptr), *tau2b(nullptr);
			tau1a = (GenParticle*) branchAll->At(a1->D1);
			tau1b = (GenParticle*) branchAll->At(a1->D2);
			tau2a = (GenParticle*) branchAll->At(a2->D1);
			tau2b = (GenParticle*) branchAll->At(a2->D2);

			// TLorentzVector tau1aMom,tau1bMom, tau2aMom, tau2bMom;
			// tau1aMom = tau1a->P4();
			// tau1bMom = tau1b->P4();
			// tau2aMom = tau2a->P4();
			// tau2bMom = tau2b->P4();

			GenParticle *charged1a = getChargedObject(branchAll, tau1a);
			GenParticle *charged1b = getChargedObject(branchAll, tau1b);
			GenParticle *charged2a = getChargedObject(branchAll, tau2a);
			GenParticle *charged2b = getChargedObject(branchAll, tau2b);
			
			// This selects events where each tau only has 1 charged product...
			// dunno what to do about evts where the tau decays into charged things *including* muon
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

					// Can use truth mu/track pairs for further analysis

					// Cleanup - Dirty!
					if (!muTruth1) delete muTruth1;
					if (!trackTruth1) delete trackTruth1;
					if (!muTruth2) delete muTruth2;
					if (!trackTruth2) delete trackTruth2;
				}
			} else {
				throw runtime_error("Not all prongs found!");
			} // end if(charged1a....)
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

		Track *candTk(nullptr);
		bool atLeastTk2p5 = false; // to monitor if theres a tk with pT > 2.5
		bool atLeastTk2p5OS = false; // same but for OS tk-muon
		for(int a = 0; a < branchTracks->GetEntries(); a++){
			candTk = (Track*) branchTracks->At(a);

			if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
				&& (candTk->PT != mu2->PT)
				&& checkTrackLoose(candTk)
			){
				// Store track in suitable vector
				fillTrackVectors(candTk, mu1, mu2, &tk1_1, &tk2_1);

				if (checkTrackTight(candTk)){
					atLeastTk2p5 = true;

					fillTrackVectors(candTk, mu1, mu2, &tk1_2p5, &tk2_2p5);

					if (checkTkMuOS(candTk, mu1)) {
						// SIGNAL REGION TRACKS
						fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS, &tk2_2p5_OS);
					}					
				}  else {
					fillTrackVectors(candTk, mu1, mu2, &tk1_1to2p5, &tk2_1to2p5);
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

	// Do any normalising, hist manipulation here...

	////////////////////////////
	// Setup filenames, paths //
	////////////////////////////
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
	
	//////////////////////////////////////////////////////////////////////
	// Plot histograms to PDF - can do anyhting that inherits from TH1 //
	//////////////////////////////////////////////////////////////////////
	// drawHistAndSave(histmu1PT, "HISTE", "muLeadingPT", directory, app);
	// drawHistAndSave(histDEtaVsDPhiMuMu, "COLZ", "DEtaVsDPhiMuMu", directory, app);

	/////////////////////////
	// Save hists to file //
	/////////////////////////
	// TFile* outFile = TFile::Open((directory+"/output_"+delph+"_"+app+".root").c_str(),"UPDATE");

	// histNMu->Write("",TObject::kOverwrite);
	// if (doSignal){
	// 	histPID->Write("",TObject::kOverwrite);
	// }

	// outFile->Close();

}
