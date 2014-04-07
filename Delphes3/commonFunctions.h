#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

// ROOT headers
#include "TStyle.h"

// BOOST headers
// Need to add 
// -I $(HOME)/boost_1_55_0 -I $(HOME)/boost_1_55_0_install/include 
// to CXXFLAGS in Delphes/Makefile
#include <boost/lexical_cast.hpp>
// #include <boost/filesystem.hpp>

using std::cout;
using std::endl;

/**
 * This header contains common functiosn for all my Delphes analysis scripts
 *
 * Robin Aggleton 2014
 */


/**
 * For splitting strings based on a delimiter
 * @param s     string to be split
 * @param delim character to be used as delimiter, e.g. ':'
 * @param elems blank vector to be used for split results
 * @return      Returns the input vector elems, filled with s split based on delim
 */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

// For sorting track vectors by pT
// Ideally we'd use the templated methods in classes/SortableObject.h ...
bool sortTracksByPT(Track* a, Track* b){ 
	return (a->PT) > (b->PT); 
}

/**
 * From the daughters of 2 taus from a given "a" Higgs, decide which is track and which is mu
 * @param  mu muon GenParticle, will be assigned to an object (a or b)
 * @param  tk track GenParticle, will be assigned to an object (a or b)
 * @param  a  One of the two gen particles from a pair of taus
 * @param  b  The other of the two GenParticles from th pair of taus
 * @return    Returns TRUE if assignemnt OK - if neither a nor b is a muon then returns FALSE
 */
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

/**
 * Get the 3 immediate duaghters of the tau (tau neutrino + qqbar or lepton+neutrino)
 * @param branchAll  The branch of ALL gen particles
 * @param tau        Tau to get decay product from
 * @return           Vector of the tau's 3 immediate daughters
 */
std::vector<GenParticle*> getTauDaughters(TClonesArray *branchAll, GenParticle *tau) { 

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

/**
 * From tau, get the final, stable, charged decay product
 * @param  branchAll The branch of ALL gen particles, to loop through decay chain
 * @param  tau       Tau to get decay product form
 * @return           Returns charged, stable, decay product of param tau
 */
GenParticle* getChargedObject(TClonesArray* branchAll, GenParticle* tau) { 

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

/**
 * Draw a histogram and save it to PDF, saves file as directory/filename_<delphes setup from directory>_app.pdf
 * @param h         Hist to draw (TObject* to handle THStacks, which don't inherit from TH1)
 * @param drawOpt   Options for drawing hist
 * @param filename  Main filename (e.g. Mu1Pt)
 * @param directory Directory for output PDF
 * @param app       Appendage eg muRand, sig
 */
void drawHistAndSave(TObject* h, std::string drawOpt, std::string filename, std::string directory, std::string app){

	TH1::SetDefaultSumw2();

	TCanvas c;
	h->Draw(drawOpt.c_str());

	std::vector<std::string> elems2;
	split(directory, '_', elems2);
	std::string delph = elems2[elems2.size()-1];
	
	c.SaveAs((directory+"/"+filename+"_"+delph+"_"+app+".pdf").c_str());
}

/**
 * Rescale histogram so integral = unity
 * @param h Poitner to histogram to be rescaled
 */
void normaliseHist(TH1* h){
	if (h->Integral() != 0){
		h->Scale(1./h->Integral());
	}
}


/**
 * Draws mass correlation plot, saves file as directory/filename_<delphes setup from directory>_app.pdf
 * @param title         Title at top of plot
 * @param histM1_0to1   Pointer to TH1 for m2 bin 0 to 1
 * @param histM1_1to2   Pointer to TH1 for m2 bin 1 to 2
 * @param histM1_2to3   Pointer to TH1 for m2 bin 2 to 3
 * @param histM1_3toInf Pointer to TH1 for m2 bin 3 to Inf
 * @param filename      Output filename
 * @param directory     Directory to put output file in
 * @param app           Any appendage onto the filename e.g. muRand
 */
void drawMassPlot(std::string title, TH1* histM1_0to1, TH1* histM1_1to2, TH1* histM1_2to3, TH1* histM1_3toInf, std::string filename, std::string directory, std::string app){

	gStyle->SetOptStat("");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(kWhite);

	TH1::SetDefaultSumw2();

	TCanvas c;
	c.SetCanvasSize(500, 600);
	TPad* pad1 = new TPad("pad1","",0,0.30,1,1);
	pad1->SetBottomMargin(0.02);
	pad1->SetRightMargin(0.05); // The ratio plot below inherits the right and left margins settings here!
	pad1->SetLeftMargin(0.15);
	pad1->Draw();
	pad1->cd();
	
	THStack histM1_M2("hM1_M2",title.c_str());
	histM1_0to1->SetLineColor(kBlack);
	histM1_0to1->SetLineWidth(2);
	normaliseHist(histM1_0to1);
	histM1_M2.Add(histM1_0to1);
	
	histM1_1to2->SetLineColor(kRed);
	histM1_1to2->SetLineWidth(2);
	normaliseHist(histM1_1to2);
	histM1_M2.Add(histM1_1to2);

	histM1_2to3->SetLineColor(kGreen+1);
	histM1_2to3->SetLineWidth(2);
	normaliseHist(histM1_2to3);
	histM1_M2.Add(histM1_2to3);

	histM1_3toInf->SetLineColor(kBlue);
	histM1_3toInf->SetLineWidth(2);
	normaliseHist(histM1_3toInf);
	histM1_M2.Add(histM1_3toInf);

	histM1_M2.Draw("nostack,HISTE");

	// Turn off x-axis label, set y-axis sizes
	// This is THE ONLY WAY
	// SetXtitle doesn't work, neither on the Stack or the individual hists
	histM1_M2.GetXaxis()->SetTitleOffset(999); // turn off x-axis title
	histM1_M2.GetXaxis()->SetLabelOffset(999); // turn off x-axis labels
	// histM1_M2.GetYaxis()->SetTitleOffset(0.08);
	histM1_M2.GetYaxis()->SetTitleSize(0.05);
	histM1_M2.Draw("nostack,HISTE");

	TLegend leg(0.75,0.67,0.93,0.88);
	leg.SetFillColor(kWhite);
	leg.AddEntry(histM1_0to1,"m_{2} = 0-1 GeV","l");
	leg.AddEntry(histM1_1to2,"m_{2} = 1-2 GeV","l");
	leg.AddEntry(histM1_2to3,"m_{2} = 2-3 GeV","l");
	leg.AddEntry(histM1_3toInf,"m_{2} > 3 GeV","l");
	leg.Draw();

	// Essential to ensure that pad2 isn't drawn INSIDE pad1	
	c.cd();

	TPad* pad2 = new TPad("pad2","",0,0,1,0.30);
	pad2->SetBottomMargin(0.25);
	pad2->SetLeftMargin(pad1->GetLeftMargin());
	pad2->SetRightMargin(pad1->GetRightMargin());
	pad2->SetTopMargin(0.05);
	pad2->SetTicks(1,1); // Puts tick marks on upper x axis and right y axis (yes, fn of Pad, not the Hist...)
	pad2->Draw();
	pad2->cd();

	TH1D* histM1_1to2_copy   = (TH1D*)histM1_1to2->Clone();
	TH1D* histM1_2to3_copy   = (TH1D*)histM1_2to3->Clone();
	TH1D* histM1_3toInf_copy = (TH1D*)histM1_3toInf->Clone();

	histM1_1to2_copy->Divide(histM1_0to1);
	histM1_2to3_copy->Divide(histM1_0to1);
	histM1_3toInf_copy->Divide(histM1_0to1);

	histM1_1to2_copy->SetTitle("");
	histM1_1to2_copy->SetXTitle("m(#mu_{1}-tk) [GeV]");
	histM1_1to2_copy->SetYTitle("#frac{i^{th} m_{2} bin}{1^{st} m_{2} bin}");
	histM1_1to2_copy->GetXaxis()->SetTitleSize(0.1);
	histM1_1to2_copy->GetXaxis()->SetTitleOffset(0.9);
	histM1_1to2_copy->GetXaxis()->SetLabelSize(0.08);

	histM1_1to2_copy->GetYaxis()->CenterTitle(kTRUE);
	histM1_1to2_copy->GetYaxis()->SetTitleSize(0.1);
	histM1_1to2_copy->GetYaxis()->SetTitleOffset(0.6);
	// histM1_1to2_copy->GetYaxis()->SetNdivisions(3,5,0);
	histM1_1to2_copy->GetYaxis()->SetLabelSize(0.08);

	histM1_1to2_copy->Draw("ep");
	histM1_2to3_copy->Draw("epSAME");
	histM1_3toInf_copy->Draw("epSAME");

	TLine *line = new TLine(0,1,10,1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw();

	// To get the delphes config used, which is always the last part of the directory name
	std::vector<std::string> elems2;
	split(directory, '_', elems2);
	std::string delph = elems2[elems2.size()-1];
	
	c.SaveAs((directory+"/"+filename+"_"+delph+"_"+app+".pdf").c_str());
	delete pad1;
	delete pad2;
}

/**
 * Add Delphes ntuples to TChain so we can process them in one go
 * @param chain    Pointer to TChain to add files to
 * @param doSignal Flag TRUE to do signalMC, FALSE to do QCD
 * @param doMu     Flag TRUE to use sample that force B hadrons to decay to mu 
 * @param doHLT    Flag TRUE to use signal sample that emulates HLT conditions (Mu17_Mu8)
 */
void addInputFiles(TChain* chain, bool doSignal, bool doMu, bool doHLT){
	// Create chain of root trees
	int nFiles = 0; // number of files to be added
	std::string stem = ""; // folder & file stem, expect files to be named like myFile_i.root, where i = 1 -> nFiles

	if (doSignal){
		if (doHLT){
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_10_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_1_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_2_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_3_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_4_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_5_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_6_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_7_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_8_HLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_9_HLT_bare.root");
			cout << "Doing signal with HLT cuts" << endl;
			stem = "Signal_1prong_500K_bare/signal_1prong_500K_HLT_";
			nFiles = 10;
		} else { 
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_10_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_1_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_2_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_3_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_4_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_5_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_6_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_7_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_8_NoHLT_bare.root");
			// chain->Add("Signal_1prong_500K_bare/signal_1prong_500K_9_NoHLT_bare.root");
			cout << "Doing signal without HLT cuts" << endl;
			stem = "Signal_1prong_500K_bare/signal_1prong_500K_NoHLT_";
			nFiles = 10;
		}
	} else {
		if (doMu){
			if(doHLT){
				cout << "Doing QCDb_mu with HLT cuts" << endl;
				stem = "QCDb_mu_pthatmin20_Mu17_Mu8_bare/QCDb_mu_pthatmin20_Mu17_Mu8_";
				nFiles = 300;
			} else{
				cout << "Doing QCDb_mu without HLT cuts" << endl;
				stem = "QCDb_mu_pthatmin20_bare/QCDb_mu_pthatmin20_";
				nFiles = 60;
			}
		} else {
			cout << "Doing QCDb" << endl;
			stem = "QCDb_cleanTk/QCDb_";
			nFiles = 10;
			// chain->Add("QCDb_cleanTk/QCDb_10.root");
			// chain->Add("QCDb_cleanTk/QCDb_2.root");
			// chain->Add("QCDb_cleanTk/QCDb_3.root");
			// chain->Add("QCDb_cleanTk/QCDb_4.root");
			// chain->Add("QCDb_cleanTk/QCDb_5.root");
			// chain->Add("QCDb_cleanTk/QCDb_6.root");
			// chain->Add("QCDb_cleanTk/QCDb_7.root");
			// chain->Add("QCDb_cleanTk/QCDb_8.root");
			// chain->Add("QCDb_cleanTk/QCDb_9.root");
		}
	}
	for (int i = 1; i <= nFiles; i ++){
		cout << "Adding " << stem+boost::lexical_cast<std::string>(i)+".root" << endl;
		chain->Add((stem+boost::lexical_cast<std::string>(i)+".root").c_str());
	}
}