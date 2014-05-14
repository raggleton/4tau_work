#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>

// ROOT headers
#include "TStyle.h"

// BOOST headers
// Need to add 
// -I $(HOME)/boost_1_55_0 -I $(HOME)/boost_1_55_0_install/include 
// to CXXFLAGS in Delphes/Makefile
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// My headers
#include "programOpts.h"

using std::cout;
using std::endl;

namespace fs = boost::filesystem;
namespace po = boost::program_options;

/**
 * This header contains common functions for all my Delphes analysis scripts
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

// Get directory that input file was in - put plots in there
std::string getDirectory(TFile* f){	
	std::string fullpath = f->GetDirectory("")->GetName();
	std::vector<std::string> elems;
	split(fullpath, '/', elems);
	return elems[0];
}

// Get Delphes file config used - last part of directory name
std::string getDelph(std::string directory){	
	std::vector<std::string> elems2;
	split(directory, '_', elems2);
	return elems2[elems2.size()-1];
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
 * Print outs cand's daughters info (index in branchAll array, PID).
 * Useful for debugging
 * @param branchAll  Branch of ALL GenParticles
 * @param cand       Particle whose daughters you wish to print to screen 
 */
void printParticleDaughters(TClonesArray *branchAll, GenParticle* cand){
	if (cand->D1 != -1 && cand->D2 != -1){
		for(int i = cand->D1; i <= cand->D2; i++){
			cout << "Daughter at " << i << ", PID: " << ((GenParticle*)branchAll->At(i))->PID << endl;
		}
	} else
		cout << "No daughters, stable." << endl;
}

/**
 * Get the 3 immediate daughters of the tau (tau neutrino + qqbar or lepton+neutrino)
 * @param branchAll  The branch of ALL gen particles
 * @param tau        Tau to get decay product from
 * @return           Vector of the tau's 3 immediate daughters
 */

// CURRENTLY DEPRECIATED

// std::vector<GenParticle*> getTauDaughters(TClonesArray *branchAll, GenParticle *tau) { 

// 	cout << "get tau daughters" << endl;

// 	std::vector<GenParticle*> tauDaughters;

// 	bool foundProduct = false;
// 	while(!foundProduct){
// 		cout << tauDaughters.size() << endl;
// 		cout << "D1: " << tau->D1 << endl;
// 		cout << "D2: " << tau->D2 << endl;
		
// 		printParticleDaughters(branchAll, tau);

// 		if(tau->D1 == tau->D2) { // tau -> tau
// 			cout << "D1 = D2" << endl;
// 			tau = (GenParticle*) branchAll->At(tau->D1);
// 		} else if ((tau->D2)-(tau->D1) == 1) { //tau->tau+gamma
			
// 			cout << "Tau -> tau + gamma" <<endl;
			
// 			if ( fabs(((GenParticle*) branchAll->At(tau->D1))->PID) == 22) {
			
// 				cout << "Tau is " << tau->D2 << " check " << ((GenParticle*) branchAll->At(tau->D2))->PID << endl;
// 				tau = (GenParticle*) branchAll->At(tau->D2);
			
// 			} else{
			
// 				cout << "Tau is " << tau->D1 << " check " << ((GenParticle*) branchAll->At(tau->D1))->PID << endl;
// 				tau = (GenParticle*) branchAll->At(tau->D1);
// 			}
		
// 		} else {
// 			for(int i = tau->D1; i <= tau->D2; i++){
// 				tauDaughters.push_back((GenParticle*) branchAll->At(i));
// 				// tauDaughters.push_back(i);
// 			}
// 			if (tauDaughters.size() ==3){
// 				cout << "3 products" << endl;
// 				foundProduct = true;
// 				return tauDaughters;
// 			}
// 		}
// 	}
// }

/**
 * From tau, get the final, stable, charged decay product.
 * Does this by stepping through "generations" - looking at daughters of all the particles in "current" vector
 * Eliminates any duplicates, (either in previoius or current generations),
 * and if the particle is final state, we return it.
 * Checks at the end to ensure only 1 charged final state particle.
 * 
 * @param  branchAll The branch of ALL gen particles, to loop through decay chain
 * @param  tau       Tau to get decay product from
 * @return           Returns charged, stable, decay product of param tau
 */
GenParticle* getChargedObject(TClonesArray* branchAll, GenParticle* tau) { 

	std::vector<GenParticle*> history; // hold all *unique* particles in the decay chain (stores position number of particle in event listing)
	std::vector<GenParticle*> current = { tau };
	std::vector<GenParticle*> next; // holds decay products, not nec. all unique
	
	GenParticle *prong = nullptr; // holds final charged product
	
	int nCharged = 0; // count charged particles - shouldn't have more than 1!
	bool foundOne = false;

	while (current.size()>0){ // if current > 0 we haven't gone through all the particles
		for (unsigned currentElem1 = 0; currentElem1 < current.size(); currentElem1++){
			
			// Check 1 - is this already in current?
			// Could probably do more efficiently using the unique function on std::vector
			bool alreadyDone = false;

			for (unsigned currentElem2 = 0; currentElem2 < currentElem1; currentElem2++){
				if ((current[currentElem1] == current[currentElem2]) && (currentElem1 != 0)) {
					// cout << "--Found a duplicate" << endl;
					alreadyDone = true;
				}
			}

			// Check 2 - is this already in history?
			if (!alreadyDone){
				// cout << "-not in current" << endl;
				for (unsigned historyElem = 0; historyElem < history.size(); historyElem++){
					if ((current[currentElem1] == history[historyElem]) && (historyElem!=0)){
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
				if (current[currentElem1]->Status == 1 || (current[currentElem1]->D1 == -1 && current[currentElem1]->D2 == -1) ) {

					// Check if charged
					if (current[currentElem1]->Charge != 0) {
						// cout << "FINAL PRONG " << current[currentElem1]->PID << endl;
						prong = current[currentElem1];
						foundOne = true;
						nCharged++;
					}

					history.push_back(current[currentElem1]);
					// cout << "pushed into history" << endl;
				} else {
					// Load its daughters no. into next
					for (int daughter = current[currentElem1]->D1; daughter <= current[currentElem1]->D2; daughter++) {
						// cout << "pushing" << endl;
						next.push_back((GenParticle*) branchAll->At(daughter));
					}

					// Load it into history
					history.push_back(current[currentElem1]);
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

	if (!foundOne) {
		cout << "No prongs!!!" << endl;
	}

	if (nCharged == 1){
		// cout << "PRONG: " << prong->PID << endl;
		return prong;
	} else {
		cout << "Damn, more than 1 prong!!!! Got " << nCharged << endl;
		return nullptr;
	}

}

/**
 * Draw a histogram (1D or 2D!) and save it to PDF, as <directory>/<filename>_<delphes setup from <directory>>_<app>.pdf
 * @param h         Hist to draw (TObject* to handle THStacks, which don't inherit from TH1)
 * @param drawOpt   Options for drawing hist
 * @param filename  Main filename (e.g. Mu1Pt)
 * @param directory Directory for output PDF
 * @param app       Appendage eg muRand, sig
 */
void drawHistAndSave(TObject* h, std::string drawOpt, std::string filename, std::string directory, std::string app, bool drawLogY = false){
	gStyle->SetOptStat("ne"); // display name and # entries only
	gStyle->SetPaintTextFormat(".3g"); // set text format to be printed
	
	TH1::SetDefaultSumw2();

	TCanvas c;
	if (drawLogY) c.SetLogy();
	h->Draw(drawOpt.c_str());

	std::string delph = getDelph(directory);
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
	histM1_1to2_copy->SetXTitle(histM1_1to2->GetXaxis()->GetTitle());
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

	// Draw dotted line at 1 on ratio plot
	double min = histM1_0to1->GetBinLowEdge(1);
	double max = histM1_0to1->GetBinLowEdge(histM1_0to1->GetNbinsX()+1);
	TLine *line = new TLine(min,1,max,1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw();

	// To get the delphes config used, which is always the last part of the directory name
	std::string delph = getDelph(directory);
	
	c.SaveAs((directory+"/"+filename+"_"+delph+"_"+app+".pdf").c_str());
	delete pad1;
	delete pad2;
}
