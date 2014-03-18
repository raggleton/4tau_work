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
// using std::cout;
// using std::endl;
using namespace std;
// For sorting track vectors by pT
// Ideally we'd use the templated methods in classes/SortableObject.h ...
// bool sortTracksByPT(Track* a, Track* b){ 
// 	return (a->PT) > (b->PT); 
// }

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

	bool doSignal = true;
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

	// Loop over all events
	Long64_t numberOfEntries = treeReader->GetEntries();
	cout << "Nevts : " << numberOfEntries <<endl;
	bool stop = false;

	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		if (branchGenMuons->GetEntries() < 2) continue; // skip if <2 muons!

		//////////////////////////////////////////////////////////////////////
		// First, get the two highest pT muons in the event, store their pT //
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

		// Vectors of tracks with pT > 1, withindR < 0.5 of respective muons + other cuts
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
		// std::sort(tk1.begin(), tk1.end(), sortTracksByPT);
		// std::sort(tk2.begin(), tk2.end(), sortTracksByPT);

		// for (unsigned a = 0; a<tk1.size(); a++) {
		// 	cout << tk1.at(a)->PT << endl;
		// }

	} // end of event loop


	// TCanvas c;
	// std::string name("");
	// std::string app("");
	// if (doSignal) {
	// 	// name = "Signal_";
	// 	name = "Signal_1prong_";
	// 	// name = "Signal_3prong_";
	// 	app = "_sig";
	// } else {
	// 	app = "_bg";
	// 	if (doMu)
	// 		name = "QCDb_mu_";
	// 	else
	// 		name = "QCDb_";
	// }
	// if (swapMuRandomly)
	// 	app += "_muRand";
	
	// app += "_samePtEta";

	// histNMu->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NMu_clean"+app+".pdf").c_str());
	// histNMu1->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NMu1_clean"+app+".pdf").c_str());
	// histNMu2->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NMu2_clean"+app+".pdf").c_str());

	// histMu1Pt->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/Mu1Pt_clean"+app+".pdf").c_str());
	// histMu2Pt->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/Mu2Pt_clean"+app+".pdf").c_str());

	// histNuPt->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NuPt_clean"+app+".pdf").c_str());
	
	// histMu1PtSel->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/Mu1PtSel_clean"+app+".pdf").c_str());
	// histMu2PtSel->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/Mu2PtSel_clean"+app+".pdf").c_str());

	// histNTracks1->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NTracks1_NS_clean"+app+".pdf").c_str());
	// histNTracks2->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NTracks2_NS_clean"+app+".pdf").c_str());

	// histNTracks1OS->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NTracks1_OS_clean"+app+".pdf").c_str());
	// histNTracks2OS->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NTracks2_OS_clean"+app+".pdf").c_str());

	// // histNTracksCum1->Draw("HISTE");
	// // c.SaveAs((name+"cleanTk/NTracks1Cum_NS_clean"+app+".pdf").c_str());
	// // histNTracksCum2->Draw("HISTE");
	// // c.SaveAs((name+"cleanTk/NTracks2Cum_NS_clean"+app+".pdf").c_str());

	// // histNTracksCum1OS->Draw("HISTE");
	// // c.SaveAs((name+"cleanTk/NTracks1Cum_OS_clean"+app+".pdf").c_str());
	// // histNTracksCum2OS->Draw("HISTE");
	// // c.SaveAs((name+"cleanTk/NTracks2Cum_OS_clean"+app+".pdf").c_str());

	// histDRMuMu->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/DRMuMu_clean"+app+".pdf").c_str());
	// histDEtaVsDPhiMuMu->Draw("COLZ");
	// c.SaveAs((name+"cleanTk/DEtaVsDPhiMuMu_clean"+app+".pdf").c_str());

	// histNTk->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NTk_clean"+app+".pdf").c_str());
	// histNTk1->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NTk1_clean"+app+".pdf").c_str());
	// histNTk25->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/NTk25_clean"+app+".pdf").c_str());

	// if (doSignal){
	// 	histDRa1->Draw("HISTE");
	// 	c.SaveAs((name+"cleanTk/DRa1_clean"+app+".pdf").c_str());
	// 	histDRa2->Draw("HISTE");
	// 	c.SaveAs((name+"cleanTk/DRa2_clean"+app+".pdf").c_str());
	// 	histPID->Draw("HISTE");
	// 	c.SaveAs((name+"cleanTk/PID_clean"+app+".pdf").c_str());
	// }

	// histTroublePt->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleTkPt_clean"+app+".pdf").c_str());
	// histTroublePID->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleTkPID_clean"+app+".pdf").c_str());
	// histTroubleEta->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleTkEta_clean"+app+".pdf").c_str());
	// histTroublePhi->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleTkPhi_clean"+app+".pdf").c_str());
	// histTroubleMatch->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleMatch_clean"+app+".pdf").c_str());
	// histTroubleDRMuMu->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleDRMuMu_clean"+app+".pdf").c_str());

	// histTroubleDPhiMuMu->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleDPhiMuMu_clean"+app+".pdf").c_str());
	// histTroubleDEtaMuMu->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleDEtaMuMu_clean"+app+".pdf").c_str());
	// histTroubleMu1Pt->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleMu1Pt_clean"+app+".pdf").c_str());
	// histTroubleMu2Pt->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/TroubleMu2Pt_clean"+app+".pdf").c_str());

	// TArc problemRing1(0,0,1,0,90);
	// problemRing1.SetLineColor(kRed);
	// problemRing1.SetLineWidth(2);
	// problemRing1.SetFillStyle(0);

	// TArc problemRing2(0,0,2,0,90);
	// problemRing2.SetLineColor(kRed);
	// problemRing2.SetLineWidth(2);
	// problemRing2.SetFillStyle(0);

	// histTroubleEtaVsPhi1->Draw("COLZ");
	// problemRing1.Draw("only");
	// c.SaveAs((name+"cleanTk/TroubleEtaVsPhi1_clean"+app+".pdf").c_str());
	// histTroubleEtaVsPhi2->Draw("COLZ");
	// problemRing1.Draw("only");
	// problemRing2.Draw("only");
	// c.SaveAs((name+"cleanTk/TroubleEtaVsPhi2_clean"+app+".pdf").c_str());

	// histRand->Draw("HISTE");
	// c.SaveAs((name+"cleanTk/RandTest"+app+".pdf").c_str());

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
