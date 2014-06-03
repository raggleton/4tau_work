#include <vector>
#include <iostream>
#include <sstream>

#include "TLegend.h"
#include "TCanvas.h" 
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLine.h"


/**
 * This just post-processes the output from massPlots program, to combine with Alexei's reuslts
 * and other stuff
 */

// template<typename T>

std::string intToString(int n) {
	ostringstream convert;
	convert << n;
	return convert.str();
}

TH2D* combinePlots(std::vector<TH2D*> plots, std::vector<double> scalingFactors) {
	TH2D* h = (TH2D*)plots[0]->Clone(plots[0]->GetName());
	h->Scale(scalingFactors[0]);
	for (unsigned i = 1; i < plots.size(); i++) {
		h->Add(plots[i], scalingFactors[i]);
	}
	return h;
}
TH1D* combinePlots(std::vector<TH1D*> plots, std::vector<double> scalingFactors) {
	TH1D* h = (TH1D*)plots[0]->Clone(plots[0]->GetName());
	h->Scale(scalingFactors[0]);
	std::vector<double> massBins;
	massBins.push_back(0);
	massBins.push_back(1);
	massBins.push_back(2);
	massBins.push_back(3);
	massBins.push_back(10);
	int nBinsX = massBins.size()-1;
	for (unsigned i = 1; i < plots.size(); i++) {
		TH1D* hTmp = (TH1D*)plots[i]->Rebin(nBinsX, plots[0]->GetTitle(), &massBins[0]);
		h->Add(hTmp, scalingFactors[i]);
	}
	return h;
}

// template<typename T>
void drawHistAndSave(TH1* h, 
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
	// if (!std::is_same<T, THStack>::value) {
		if (drawOpt.find("TEXT") != std::string::npos) {
			h->SetMarkerSize(1.5*h->GetMarkerSize());
		} else {
			h->SetMarkerStyle(20);
			h->SetMarkerColor(h->GetLineColor());
		}
	// }

	h->Draw(drawOpt.c_str());

	// Save to PDF
	std::string sep("");
	if (directory != "") sep = "/";
	std::string sep2("");
	if (app != "") sep2 = "_";

	c.SaveAs((directory+sep+filename+sep2+app+".pdf").c_str());
}

void massPlotPost() {
	gStyle->SetOptStat(""); // DOES NOTHING AS HIST ALREADY HAS OPT STATS!!!!
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(kWhite);

	TH1::SetDefaultSumw2();
	TCanvas* c1 = new TCanvas("c1");

	////////////////////////////
	// Open files, get hists //
	////////////////////////////
	TFile* fQCDb = TFile::Open("../QCDb_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* fQCDScatter = TFile::Open("../QCDScatter_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* fQCDc = TFile::Open("../QCDc_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");

	std::vector<TFile*> files;
	files.push_back(fQCDb);
	files.push_back(fQCDScatter);
	files.push_back(fQCDc);

	TH1D* histCorr1D_side_1to2p5_Robin_QCDb = (TH1D*) fQCDb->Get("hCorr1D_side_1to2p5");
	TH1D* histCorr1D_side_1to2p5_Robin_QCDScatter = (TH1D*) fQCDScatter->Get("hCorr1D_side_1to2p5");
	TH1D* histCorr1D_side_1to2p5_Robin_QCDc = (TH1D*) fQCDc->Get("hCorr1D_side_1to2p5");
	histCorr1D_side_1to2p5_Robin_QCDb->SetStats(kFALSE);
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetStats(kFALSE);
	histCorr1D_side_1to2p5_Robin_QCDc->SetStats(kFALSE);
	
	////////////////////////////
	// Make Alexei's 1D hist //
	////////////////////////////
	TH1D* histCorr1D_side_1to2p5_Alexei = (TH1D*) histCorr1D_side_1to2p5_Robin_QCDb->Clone("hCorr1D_side_1to2p5_Alexei");
	double arrContents[] = {1.00, 1.02, 0.95, 1.21, 0.96, 0.93, 0.97, 1.1, 1.1, 1.03};
	double arrErrors[] = {0.05, 0.03, 0.05, 0.08, 0.04, 0.05, 0.07, 0.13, 0.12, 0.22};
	for (unsigned i = 1; i <= histCorr1D_side_1to2p5_Robin_QCDb->GetNbinsX(); i++) {
		histCorr1D_side_1to2p5_Alexei->SetBinContent(i,arrContents[i-1]);
		histCorr1D_side_1to2p5_Alexei->SetBinError(i,arrErrors[i-1]);
	}

	//////////////////////////////////////////
	// Plot individual components together //
	//////////////////////////////////////////
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetLineColor(kGreen+2);
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetMarkerColor(kGreen+2);
	histCorr1D_side_1to2p5_Robin_QCDc->SetLineColor(kOrange+2);
	histCorr1D_side_1to2p5_Robin_QCDc->SetMarkerColor(kOrange+2);
	histCorr1D_side_1to2p5_Alexei->SetLineColor(kRed);
	histCorr1D_side_1to2p5_Alexei->SetMarkerColor(kRed);
	THStack stack("stack","");
	stack.Add(histCorr1D_side_1to2p5_Robin_QCDb);
	stack.Add(histCorr1D_side_1to2p5_Robin_QCDScatter);
	// stack.Add(histCorr1D_side_1to2p5_Robin_QCDc);
	stack.Add(histCorr1D_side_1to2p5_Alexei);
	stack.Draw("EPNOSTACK");
	(stack.GetHistogram())->SetXTitle("Bin");
	(stack.GetHistogram())->SetYTitle("Correlation coefficient");

	// Add a legend
	TLegend leg(0.65,0.7,0.88,0.88);
	leg.SetFillColor(kWhite);
	leg.AddEntry(histCorr1D_side_1to2p5_Robin_QCDb,"QCDb","lp");
	// leg.AddEntry(histCorr1D_side_1to2p5_Robin_QCDc,"QCDc","lp");
	leg.AddEntry(histCorr1D_side_1to2p5_Robin_QCDScatter,"QCD q-g Scatter","lp");
	leg.AddEntry(histCorr1D_side_1to2p5_Alexei,"Data (Alexei)","lp");
	leg.Draw();

	// Draw a horizontal line at 1
	double min = histCorr1D_side_1to2p5_Alexei->GetBinLowEdge(1);
	double max = histCorr1D_side_1to2p5_Alexei->GetBinLowEdge(histCorr1D_side_1to2p5_Alexei->GetNbinsX()+1);
	TLine *line = new TLine(min,1,max,1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw();

	// Default canvas has name c1
	c1->SetTicks(1,1); // Put tick marks on top x and right y axes
	c1->SaveAs("histCorr1D_side_1to2p5.pdf");

	//////////////////////////////////////////////////////////////////
	// Now combine QCD plots - need to reweight for cross-sections //
	//////////////////////////////////////////////////////////////////

	std::vector<TH2D*> plots2D;
	std::vector<TH1D*> plots1D;

	for (unsigned i = 0; i < files.size(); i++) {
		plots2D.push_back((TH2D*)files[i]->Get("hM1vsM2_side_1to2p5"));
		plots1D.push_back((TH1D*)files[i]->Get("hM_side_1to2p5"));
	}

	std::vector<double> scalingFactors;
	scalingFactors.push_back(1.); // QCDb
	scalingFactors.push_back(1.); // QCDc
	scalingFactors.push_back(1.); // QCD scatter

	// check we haven't fluffed up vectors.
	if (plots2D.size() != scalingFactors.size()) exit(-1);

	// Create combination 2D plot (numerator)
	TH2D* histM1vsM2_side_1to2p5 = (TH2D*)combinePlots(plots2D, scalingFactors);

	// Create combination 1D sideband plot (denominator)
	TH1D* histM_side_1to2p5 = (TH1D*)combinePlots(plots1D, scalingFactors);
	
	std::vector<double> massBins;
	massBins.push_back(0);
	massBins.push_back(1);
	massBins.push_back(2);
	massBins.push_back(3);
	massBins.push_back(10);

	int nBinsX = massBins.size()-1;

	// Create 2D from 1D x 1D
	TH2D* histM1timesM1_side_1to2p5 = new TH2D("hM1timesM2_side_1to2p5","m(sideband) #times m(sideband) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	for(int a = 1; a <= nBinsX; a++){
		for (int b = 1; b <=nBinsX; b++){
			histM1timesM1_side_1to2p5->SetBinContent(a,b,histM_side_1to2p5->GetBinContent(a)*histM_side_1to2p5->GetBinContent(b));
			histM1timesM1_side_1to2p5->SetBinError(a,b,sqrt(pow(histM_side_1to2p5->GetBinContent(b)*histM_side_1to2p5->GetBinError(a),2)
															+pow(histM_side_1to2p5->GetBinContent(a)*histM_side_1to2p5->GetBinError(b),2)));
		}
	}

	// Calculate new correlation coeffs & plot
	TH2D* histM1vsM2_correlations_side_1to2p5 = (TH2D*)histM1vsM2_side_1to2p5->Clone("hM1vsM2_correlations_side_1to2p5");
	histM1vsM2_correlations_side_1to2p5->SetTitle(
		"m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband), (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_side_1to2p5->Divide(histM1timesM1_side_1to2p5);

	// Make 1D plots of unique bins from 2D correlation plot
	int nUniqueBins = (nBinsX+1)*nBinsX/2.;
	TH1D* histCorr1D_side_1to2p5_combo = new TH1D("hCorr1D_side_1to2p5_combo",";Bin;Correlation coefficient",nUniqueBins,0,nUniqueBins);
	int counter = 1;
	for (int i = 1; i <= nBinsX; i++) {
		for (int j = i; j <= nBinsX; j++) {
			std::string binLabel = "(" + intToString(i) + "," + intToString(j) + ")";
			histCorr1D_side_1to2p5_combo->SetBinContent(counter,histM1vsM2_correlations_side_1to2p5->GetBinContent(i,j));
			histCorr1D_side_1to2p5_combo->SetBinError(counter,histM1vsM2_correlations_side_1to2p5->GetBinError(i,j));
			histCorr1D_side_1to2p5_combo->GetXaxis()->SetBinLabel(counter,binLabel.c_str());
			
			counter++;
		}
	}

	std::string directory = "";
	std::string app = "";
	drawHistAndSave(histM1vsM2_side_1to2p5, "colzTEXTE","M1vsM2_side_1to2p5", directory, app);
	drawHistAndSave(histM1timesM1_side_1to2p5, "colzTEXTE","M1timesM1_side_1to2p5", directory, app);
	drawHistAndSave(histM_side_1to2p5, "HISTE", "M_side_1to2p5", directory, app);
	drawHistAndSave(histM1vsM2_correlations_side_1to2p5, "colzTEXTE","M1vsM2_correlations_side_1to2p5", directory, app);
	drawHistAndSave(histCorr1D_side_1to2p5_combo, "e1", "Correlations1D_side_1to2p5", directory, app);

	// Plot alongside Alexei's
	THStack stack2("stack2","");
	stack2.Add(histCorr1D_side_1to2p5_combo);
	stack2.Add(histCorr1D_side_1to2p5_Alexei);
	stack2.Draw("EPNOSTACK");
	(stack2.GetHistogram())->SetXTitle("Bin");
	(stack2.GetHistogram())->SetYTitle("Correlation coefficient");

	// Add a legend
	TLegend leg2(0.65,0.7,0.88,0.88);
	leg2.SetFillColor(kWhite);
	leg2.AddEntry(histCorr1D_side_1to2p5_combo,"QCD b + q-g scatter","lp");
	// leg2.AddEntry(histCorr1D_side_1to2p5_combo_Robin_QCDc,"QCDc","lp");
	leg2.AddEntry(histCorr1D_side_1to2p5_Alexei,"Data (Alexei)","lp");
	leg2.Draw();
	
	c1->SetTicks(1,1);
	c1->SaveAs("histCorr1D_side_1to2p5_combo_allQCD.pdf");
}
