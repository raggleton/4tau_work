#include "commonFunctions.h"

using std::cout;
using std::endl;

/**
 * Plot IP of tracks, and displacement of object from PV
 */

template<typename T>
 double calculateIP(T* candTk){
	// Calculate impact paramter (IP) of track object.
	// Let x = position vector to object from PV
	// Let v = velocity vector of object
	// IP = |x|*sin(theta)
	// Where cos(theta) = v.x/|v||x|
	TVector3 v = (candTk->P4()).Vect();
	TVector3 x(candTk->X,candTk->Y,candTk->Z);
	double theta = TMath::ACos(v.Unit() * x.Unit());
	return fabs(x.Mag() * TMath::Sin(theta));
}

void IP(int argc, char* argv[])
{
	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");


	ProgramOpts pOpts(argc, argv);

	// MCsource source = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	bool doSignal = pOpts.getSignal(); // for signal or not
	// bool doMu = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly = pOpts.getMuOrdering(); // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT = pOpts.getHLT(); // whether to use MC that has HLT cuts already applied or not.
	
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
	// TClonesArray *branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray *branchAll      = treeReader->UseBranch("AllParticle");

	//////////////////////
	// Book histograms //
	//////////////////////
	TH1D *histIPTracks        = new TH1D("hIPTracks" ,"Track IP, p_{T}(trk)>2.5 GeV, muon selection;Track IP [mm]; A.U.", 50,0,0.5);
	TH1D *histIPTracksOS      = new TH1D("hIPTracksOS" ,"Track IP, OS to #mu, p_{T}(trk)>2.5 GeV, muon selection;Track IP [mm]; A.U.", 50,0,0.5);
	TH1D *histIPTracksTruth   = new TH1D("hIPTracksTruth" ,"Track IP for #tau decay products, MC truth, p_{T}(trk)>2.5 GeV, no muon selection;Track IP [mm]; A.U.", 50,0,0.5);
	TH1D *histIPTracksTruthOS = new TH1D("hIPTracksTruthOS" ,"Track IP for #tau decay products, MC truth, OS to #mu, p_{T}(trk)>2.5 GeV, no muon selection;Track IP [mm]; A.U.", 50,0,0.5);

	///////////////////////
	// Loop over events //
	///////////////////////
	Long64_t numberOfEntries = getNumberEvents(treeReader, &pOpts);
	cout << "Running over " << numberOfEntries << " events" << endl;

	bool stop = false; // used to stop the loop, for debugging/testing

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
					// cout << "Problem, no truth mu for 1 and/or 2!" << endl;
				} else { 
					
/*					// Do m1 distribution in bins of m2 - for MC truth (is it actually correlated?)
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
*/
					// cout << m1 << "     " << m2 << endl;
					// if(m2 < 1.)
					// 	histM1_truth_0to1->Fill(m1);
					// else if (m2 < 2.)
					// 	histM1_truth_1to2->Fill(m1);
					// else if (m2 < 3.)
					// 	histM1_truth_2to3->Fill(m1);
					// else
					// 	histM1_truth_3toInf->Fill(m1);
					
					// Do IP for tau decay products
					double IP = calculateIP(trackTruth1);
					histIPTracksTruth->Fill(IP);
					if((trackTruth1->Charge * muTruth1->Charge)< 0) histIPTracksTruthOS->Fill(IP);
					IP = calculateIP(trackTruth2);
					histIPTracksTruth->Fill(IP);
					if((trackTruth2->Charge * muTruth2->Charge)< 0) histIPTracksTruthOS->Fill(IP); 
					// note: OS plto should be identical to non OS one, if your MC is right!
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
		// so tk1 is the track nearest to muon1, (may or may not be highest pT, depends if random swapping is on)
		// std::vector<Track*> tk1_1;
		// std::vector<Track*> tk2_1;

		// same but with pT >2.5
		// std::vector<Track*> tk1_2p5;
		// std::vector<Track*> tk2_2p5;

		// same but with pT > 2.5, OS to muon
		// std::vector<Track*> tk1_2p5_OS;
		// std::vector<Track*> tk2_2p5_OS;

		Track *candTk(0);
		for(int a = 0; a < branchTracks->GetEntries(); a++){
			candTk = (Track*) branchTracks->At(a);

			if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
				&& (candTk->PT != mu2->PT)
				&& (candTk->PT > 2.5)
				&& (fabs(candTk->Eta)<3)
			){
				// Calculate impact paramter (IP) of track object.
				// Let x = position vector to object from PV
				// Let v = velocity vector of object
				// IP = |x|*sin(theta)
				// Where cos(theta) = v.x/|v||x|
				// TVector3 v = (candTk->P4()).Vect();
				// TVector3 x(candTk->X,candTk->Y,candTk->Z);
				// double theta = TMath::ACos(v.Unit() * x.Unit());
				// double IP = x.Mag() * TMath::Sin(theta);
				double IP = calculateIP(candTk);
				histIPTracks->Fill(IP);

				// Do OS tracks (SS muon requirement is above)
				if((candTk->Charge * mu1->Charge)< 0) histIPTracksOS->Fill(IP);
			}

			// Standard track conditions
/*			if (   (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
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
			} // End of track selection*/
		} // End of track loop

		// Now pT order the track collections
		// std::sort(tk1.begin(), tk1.end(), sortTracksByPT);
		// std::sort(tk2.begin(), tk2.end(), sortTracksByPT);

		// Can now deal with signal or sideband regions
		
		// SIGNAL SELECTION
		// if (tk1_1.size() == 1 && tk2_1.size() == 1 
		// && tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {
			// do something in signal region...
		// }		

		// Clean up
		// tk1_1.clear();
		// tk2_1.clear();
		// tk1_2p5.clear();
		// tk2_2p5.clear();
		// tk1_2p5_OS.clear();
		// tk2_2p5_OS.clear();
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
	
	normaliseHist(histIPTracks);
	normaliseHist(histIPTracksOS);
	normaliseHist(histIPTracksTruth);
	normaliseHist(histIPTracksTruthOS);

	drawHistAndSave(histIPTracks, "HISTE", "TrackIP", directory, app, true);
	drawHistAndSave(histIPTracksOS, "HISTE", "TrackIPOS", directory, app, true);
	drawHistAndSave(histIPTracksTruth, "HISTE", "TrackIPTruth", directory, app, true);
	drawHistAndSave(histIPTracksTruthOS, "HISTE", "TrackIPTruthOS", directory, app, true);

	TFile* outFile = TFile::Open((directory+"/output_"+delph+"_"+app+".root").c_str(),"UPDATE");

	histIPTracks->Write("",TObject::kOverwrite);
	histIPTracksOS->Write("",TObject::kOverwrite);
	if (doSignal){
		histIPTracksTruth->Write("",TObject::kOverwrite);
		histIPTracksTruthOS->Write("",TObject::kOverwrite);
	}

	outFile->Close();
}
