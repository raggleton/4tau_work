#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <type_traits>

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
// #include "tdrstyle.C"

using std::cout;
using std::endl;

namespace fs = boost::filesystem;
namespace po = boost::program_options;

/**
 * This header contains common functions for all my Delphes analysis scripts.
 *
 * For classes/functions pertaining to user options, see programOpts.h
 * For cuts on tracks and muons, see cuts.h
 * 
 * Robin Aggleton 2014
 */

/**
 * For splitting strings based on a delimiter
 * @param s     string to be split
 * @param delim character to be used as delimiter, e.g. ':'
 * @param elems blank vector to be used for split results
 * @return      Returns input vec <elems>, filled with s split based on delim
 */
std::vector<std::string> &split(const std::string &s, 
								char delim, 
								std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

// Get directory that input file was in - put plots in there
std::string getDirectory(TFile* f) {	
	std::string fullpath = f->GetDirectory("")->GetName();
	std::vector<std::string> elems;
	split(fullpath, '/', elems);
	return elems[0];
}

// Get Delphes file config used - last part of directory name
std::string getDelph(std::string directory) {	
	std::vector<std::string> elems2;
	split(directory, '_', elems2);
	return elems2[elems2.size()-1];
}

/**
 * Dermine number of events to run over, based on user input.
 * If user inputs "-1", it will return all events in tree.
 * Note that if the user doesn't specify anything, 
 * by default it returns the total number of events in tree.
 * @param  treeReader The Delphes tree, used to determine total # events
 * @param  pOpts      ProgramOpts object that contains user's options
 * @return            Returns either total numbers of events in the tree, 
 *                    or the number of events the user specified.
 */
Long64_t getNumberEvents(ExRootTreeReader *treeReader, ProgramOpts* pOpts) {
	int userNEvents = pOpts->getNEvents();
	int treeNEvents = treeReader->GetEntries();
	if (userNEvents == -1) {
		return treeNEvents;
	} else {
		if (userNEvents > treeNEvents) {
			cout << "WARNING: You wanted " << userNEvents << " but only " 
			<< treeNEvents << " in files, so doing those ones." << endl;
			return treeNEvents;
		} else {
			return userNEvents;
		}
	}
}

/**
 * Compare two object by their pT. 
 * @param  a First obj  to compare
 * @param  b Other obj to compare
 * @return   Returns TRUE if pT of obj a > pT obj b, FALSE otherwise
 */
template <typename T>
bool sortByPT(T* a, T* b) {
	if(a->PT > b->PT) {
		return true; 
	} else {
		return false;
	}
}

/**
 * From the daughters of 2 taus from a given "a" Higgs, 
 * decide which is track and which is mu
 * @param  mu muon GenParticle, will be assigned to an object (a or b)
 * @param  tk track GenParticle, will be assigned to an object (a or b)
 * @param  a  One of the two gen particles from a pair of taus
 * @param  b  The other of the two GenParticles from th pair of taus
 * @return    Returns TRUE if assignemnt OK - 
 *            if neither a nor b is a muon then returns FALSE
 */
bool assignMuonAndTrack(GenParticle* &mu, 
						GenParticle* &tk, 
						GenParticle &a, 
						GenParticle &b) {
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
	if (aIsMu && bIsMu) {
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
void printParticleDaughters(TClonesArray *branchAll, GenParticle* cand) {
	if (cand->D1 != -1 && cand->D2 != -1) {
		for(int i = cand->D1; i <= cand->D2; i++) {
			cout << "Daughter at " << i << ", PID: " << 
			((GenParticle*)branchAll->At(i))->PID << endl;
		}
	} else
		cout << "No daughters, stable." << endl;
}

/**
 * From tau, get the final, stable, charged decay product.
 * Does this by stepping through "generations" 
 * - looking at daughters of all the particles in "curr" vector
 * Eliminates any duplicates, (either in previoius or current generations),
 * and if the particle is final state, we return it.
 * Checks at the end to ensure only 1 charged final state particle.
 * 
 * @param  branchAll The branch of ALL gen particles, to go through decay chain
 * @param  tau       Tau to get decay product from
 * @return           Returns charged, stable, decay product of param tau
 */
GenParticle* getChargedObject(TClonesArray* branchAll, GenParticle* tau) { 
	
	// hold all *unique* particles in the decay chain 
	// (stores position number of particle in event listing)
	std::vector<GenParticle*> history; 
	
	std::vector<GenParticle*> curr = { tau };
	
	// holds decay products from particles in curr, not necessarily all unique
	std::vector<GenParticle*> next; 
	
	GenParticle *prong = nullptr; // holds final charged product
	
	int nCharged = 0; // count charged particles - shouldn't have more than 1!
	bool foundOne = false;

	while (curr.size()>0) { // if 0 we haven't gone through all the particles
		// for (unsigned currElem1 = 0; currElem1 < curr.size(); currElem1++) {
		for (unsigned currElem1 = 0; currElem1 < curr.size(); currElem1++) {
			
			// Check 1 - is this already in curr?
			// Could probably do more efficiently using the unique function
			bool alreadyDone = false;

			for (unsigned currElem2 = 0; currElem2 < currElem1; currElem2++) {
				if ((curr[currElem1] == curr[currElem2]) && (currElem1 != 0)) {
					alreadyDone = true;
				}
			}

			// Check 2 - is this already in history?
			if (!alreadyDone) {
				// for (unsigned historyElem = 0; historyElem < history.size(); historyElem++) {
				for (auto historyElem : history) {
					if ((curr[currElem1] == historyElem) && (historyElem!=history[0])) {
						alreadyDone = true;
					}
				}
			}

			// Check 3 - is this final state?
			if (!alreadyDone) {
				// Check if final state. Either status 1, or D1 == D2 == -1
				if (curr[currElem1]->Status == 1 || 
					(curr[currElem1]->D1 == -1 && curr[currElem1]->D2 == -1)) {

					// Check if charged
					if (curr[currElem1]->Charge != 0) {
						prong = curr[currElem1];
						foundOne = true;
						nCharged++;
					}

					history.push_back(curr[currElem1]);
				} else {
					// Load its daughters no. into next
					for (int daughter = curr[currElem1]->D1; 
						 daughter <= curr[currElem1]->D2; daughter++) {
						next.push_back((GenParticle*) branchAll->At(daughter));
					}

					// Load it into history
					history.push_back(curr[currElem1]);
				}

			}// end of alreadyDone
		} // end of loop over curr

		// Clear curr - don't need to do as = assignment auto does this
		// curr.clear();
		// Copy next into curr
		curr = next;
		// Empty next
		next.clear();
	} // end of while(!done)

	// Just in case
	history.clear();
	curr.clear();
	next.clear();

	if (!foundOne) {
		cout << "No prongs!!!" << endl;
	}

	if (nCharged == 1) {
		return prong;
	} else {
		cout << "Damn, more than 1 prong!!!! Got " << nCharged << endl;
		return nullptr;
	}

}

/**
 * Draw a histogram (1D or 2D!) and save it to PDF, 
 * as <directory>/<filename>_<delphes setup from <directory>>_<app>.pdf
 * @param h         Hist to draw. Since THStack and TH1 don't share 
 *                  SetMarkerSize method, need to check object type. 
 *                  Suppose I should use specialisation here instead...
 * @param drawOpt   Options for drawing hist
 * @param filename  Main filename (e.g. Mu1Pt)
 * @param directory Directory for output PDF
 * @param app       Appendage eg muRand, sig
 */
template<typename T>
void drawHistAndSave(T* h, 
					 std::string drawOpt, 
					 std::string filename, 
					 std::string directory, 
					 std::string app, 
					 bool drawLogY = false) {

	gStyle->SetOptStat(""); // display name and # entries only
	gStyle->SetPaintTextFormat(".3g"); // set text format to be printed
	
	// setTDRStyle();

	TH1::SetDefaultSumw2();

	TCanvas c;
	if (drawLogY) c.SetLogy();
	
	// For plots drawing text values of bins, make text bigger
	if (!std::is_same<T, THStack>::value) {
		if (drawOpt.find("TEXT") != std::string::npos) {
			h->SetMarkerSize(1.5*h->GetMarkerSize());
		} else {
			h->SetMarkerStyle(20);
			h->SetMarkerColor(h->GetLineColor());
		}
	}

	h->Draw(drawOpt.c_str());
	if (std::is_same<T, TH1D>::value) {
		h->SetMinimum(0);
	}

	// To get the delphes config used, 
	// which is always the last part of the directory name
	std::string delph = getDelph(directory);
	
	// Save to PDF
	c.SaveAs((directory+"/"+filename+"_"+delph+"_"+app+".pdf").c_str());
}

/**
 * Rescale histogram so integral = unity
 * @param h Pointer to histogram to be rescaled
 */
void normaliseHist(TH1* h) {
	if (h->Integral() != 0) {
		h->Scale(1./h->Integral());
	}
}

/**
 * Print out integral of TH1 object (or dervied class)
 * @param h Print integral of this hist
 */
void printIntegral(TH1* h) {
	cout << h->GetName() << " : [" << h->GetTitle() << "] = " << h->Integral() << endl;
}

/**
 * Draws mass correlation plot, saves file as 
 * directory/filename_<delphes setup from directory>_app.pdf
 * @param title         Title at top of plot
 * @param histM1_0to1   Pointer to TH1 for m2 bin 0 to 1
 * @param histM1_1to2   Pointer to TH1 for m2 bin 1 to 2
 * @param histM1_2to3   Pointer to TH1 for m2 bin 2 to 3
 * @param histM1_3toInf Pointer to TH1 for m2 bin 3 to Inf
 * @param filename      Output filename
 * @param directory     Directory to put output file in
 * @param app           Any appendage onto the filename e.g. muRand
 */
void drawMassPlot(std::string title, 
				  TH1* histM1_0to1, 
				  TH1* histM1_1to2, 
				  TH1* histM1_2to3, 
				  TH1* histM1_3toInf, 
				  std::string filename, 
				  std::string directory, 
				  std::string app) {

	gStyle->SetOptStat("");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(kWhite);
	// setTDRStyle();

	TH1::SetDefaultSumw2();

	TCanvas c;
	c.SetCanvasSize(500, 600);
	TPad* pad1 = new TPad("pad1","",0,0.30,1,1);
	pad1->SetBottomMargin(0.02);
	// The ratio plot below inherits the right and left margins settings here!
	pad1->SetRightMargin(0.05);
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
	// Now do the ratio plot below
	TPad* pad2 = new TPad("pad2","",0,0,1,0.30);
	pad2->SetBottomMargin(0.25);
	pad2->SetLeftMargin(pad1->GetLeftMargin());
	pad2->SetRightMargin(pad1->GetRightMargin());
	pad2->SetTopMargin(0.05);
	pad2->SetTicks(1,1); // Put tick marks on upper x and right y axes (fn of Pad, not the Hist!)
	pad2->Draw();
	pad2->cd();

	TH1D* histM1_1to2_copy = (TH1D*)histM1_1to2->Clone();
	TH1D* histM1_2to3_copy = (TH1D*)histM1_2to3->Clone();
	TH1D* histM1_3toInf_copy = (TH1D*)histM1_3toInf->Clone();

	histM1_1to2_copy->Divide(histM1_0to1);
	histM1_2to3_copy->Divide(histM1_0to1);
	histM1_3toInf_copy->Divide(histM1_0to1);

	THStack hist_ratios("hRatios","");
	hist_ratios.Add(histM1_1to2_copy);
	hist_ratios.Add(histM1_2to3_copy);
	hist_ratios.Add(histM1_3toInf_copy);
	hist_ratios.Draw("epNOSTACK"); // Need to draw *before* using GetHistogram()

	(hist_ratios.GetHistogram())->SetXTitle(histM1_1to2->GetXaxis()->GetTitle());
	(hist_ratios.GetHistogram())->GetXaxis()->SetTitleSize(0.1);
	(hist_ratios.GetHistogram())->GetXaxis()->SetTitleOffset(0.9);
	(hist_ratios.GetHistogram())->GetXaxis()->SetLabelSize(0.08);

	(hist_ratios.GetHistogram())->SetYTitle("#frac{i^{th} m_{2} bin}{1^{st} m_{2} bin}");
	(hist_ratios.GetHistogram())->GetYaxis()->CenterTitle(kTRUE);
	(hist_ratios.GetHistogram())->GetYaxis()->SetTitleSize(0.1);
	(hist_ratios.GetHistogram())->GetYaxis()->SetTitleOffset(0.6);
	(hist_ratios.GetHistogram())->GetYaxis()->SetLabelSize(0.08);
	hist_ratios.Draw("epNOSTACK");

	// Draw dotted line at 1 on ratio plot
	double min = histM1_0to1->GetBinLowEdge(1);
	double max = histM1_0to1->GetBinLowEdge(histM1_0to1->GetNbinsX()+1);
	TLine *line = new TLine(min,1,max,1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw();

	// To get the delphes config used, 
	// which is always the last part of the directory name
	std::string delph = getDelph(directory);
	
	// Save to PDF
	c.SaveAs((directory+"/"+filename+"_"+delph+"_"+app+".pdf").c_str());
	
	delete pad1;
	delete pad2;
}

#endif