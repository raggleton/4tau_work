#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <cmath>

#include "TLegend.h"
#include "TCanvas.h" 
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPaveText.h"

using std::cout;
using std::endl;

/**
 * This just post-processes the output from massPlots program, to combine with Data reuslts
 * and other stuff
 */

const std::vector<double> massBins {0,1,2,3,10};
const int nBinsX = massBins.size()-1;

// template<typename T>
/**
 * Rescale histogram so integral = unity
 * @param h Pointer to histogram to be rescaled
 */
void normaliseHist(TH1* h) {
	TH1::SetDefaultSumw2();

	if (h->Integral() != 0) {
		h->Scale(1./h->Integral());
	}
}

std::string intToString(int n) {
	std::ostringstream convert;
	convert << n;
	return convert.str();
}

TH2D* combinePlots(std::vector<TH2D*> plots, std::vector<double> scalingFactors) {
	TH1::SetDefaultSumw2();

	TH2D* h = (TH2D*)plots[0]->Clone(plots[0]->GetName());
	h->Scale(scalingFactors[0]);
	for (unsigned i = 1; i < plots.size(); i++) {
		h->Add(plots[i], scalingFactors[i]);
	}
	return h;
}
TH1D* combinePlots(std::vector<TH1D*> plots, std::vector<double> scalingFactors) {
	TH1::SetDefaultSumw2();

	TH1D* h = (TH1D*)plots[0]->Clone(plots[0]->GetName());
	h->Scale(scalingFactors[0]);
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
					 bool drawLogY = false,
					 double ymax = -1,
					 double ymin = -1) {

	gStyle->SetOptStat(""); // display name and # entries only
	gStyle->SetPaintTextFormat(".3g"); // set text format to be printed
	gStyle->SetHistLineWidth(2);
	// setTDRStyle();

	TH1::SetDefaultSumw2();

	TCanvas c;
	if (drawLogY) c.SetLogy();
	
	// For plots drawing text values of bins, make text bigger
	// if (!std::is_same<T, THStack>::value) {
		if (drawOpt.find("TEXT") != std::string::npos) {
			// h->SetMarkerSize(1.5*h->GetMarkerSize());
		} else {
			h->SetMarkerStyle(20);
			h->SetMarkerColor(h->GetLineColor());
		}
	// }
	if (ymin != -1)
		h->SetMinimum(ymin);
	else
		h->SetMinimum(0);
	if (ymax != -1)
		h->SetMaximum(ymax);

	h->Draw(drawOpt.c_str());

	// Save to PDF
	std::string sep("");
	if (directory != "") sep = "/";
	std::string sep2("");
	if (app != "") sep2 = "_";

	c.SaveAs((directory+sep+filename+sep2+app+".pdf").c_str());
}

TH2D* create2Dfrom1D(std::vector<double> bins, TH1D* h) {
	int nBins = bins.size()-1;
	TH2D* hist2D = new TH2D("", "", nBins,&bins[0],nBins,&bins[0]);
	for(int a = 1; a <= nBins; a++){
		for (int b = 1; b <=nBins; b++){
			hist2D->SetBinContent(a,b,h->GetBinContent(a)*h->GetBinContent(b));
			hist2D->SetBinError(a,b,sqrt(pow(h->GetBinContent(b)*h->GetBinError(a),2)
															+pow(h->GetBinContent(a)*h->GetBinError(b),2)));
		}
	}
	return hist2D;
}

TH1D* unique1DBinsFrom2D(TH2D* h2D, int nUniqueBins) {
	TH1D* h1D = new TH1D("","",nUniqueBins,0,nUniqueBins);
	int counter = 1;
	for (int i = 1; i <= nBinsX; i++) {
		for (int j = i; j <= nBinsX; j++) {
			std::string binLabel = "(" + intToString(i) + "," + intToString(j) + ")";
			h1D->SetBinContent(counter,h2D->GetBinContent(i,j));
			h1D->SetBinError(counter,h2D->GetBinError(i,j));
			h1D->GetXaxis()->SetBinLabel(counter,binLabel.c_str());
			counter++;
		}
	}
	return h1D;
}

void setupCorrStuff(TH1* h){
	h->SetXTitle("Bin");
	h->SetTitleSize(0.05,"X");
	h->SetLabelSize(0.07,"X");
	h->SetLabelSize(0.05,"Y");
	h->SetYTitle("Correlation coefficient");
	h->SetTitleSize(0.05,"Y");
	h->SetMaximum(2.0);
	h->SetMinimum(0);
}

// int massPlotPost() {
int main() {
	gStyle->SetOptStat(""); // DOES NOTHING AS HIST ALREADY HAS OPT STATS!!!!
	gStyle->SetLegendBorderSize(0);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLegendFillColor(kWhite);

	TH1::SetDefaultSumw2();
	TCanvas* c1 = new TCanvas("c1");
	std::string directory = "../Combined";
	std::string app = "";

	////////////////////////////
	// Open files, get hists //
	////////////////////////////
	TFile* fQCDb = TFile::Open("../QCDb_HLT_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* fQCDScatter = TFile::Open("../QCDbcScatter_HLT_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* fSignal = TFile::Open("../Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT_dR1.root","READ");

	std::vector<TFile*> files;
	files.push_back(fQCDb);
	files.push_back(fQCDScatter);

	TH1D* histCorr1D_side_1to2p5_Robin_QCDb = (TH1D*) fQCDb->Get("hCorr1D_side_1to2p5");
	TH1D* histCorr1D_side_1to2p5_Robin_QCDScatter = (TH1D*) fQCDScatter->Get("hCorr1D_side_1to2p5");
	TH1D* histCorr1D_side_1to2p5_Robin_Signal = (TH1D*) fSignal->Get("hCorr1D_side_1to2p5");
	histCorr1D_side_1to2p5_Robin_QCDb->SetStats(kFALSE);
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetStats(kFALSE);
	histCorr1D_side_1to2p5_Robin_Signal->SetStats(kFALSE);
	
	////////////////////////////
	// Make Data 1D hist //
	////////////////////////////
	TH1D* histCorr1D_side_1to2p5_Data = (TH1D*) histCorr1D_side_1to2p5_Robin_QCDb->Clone("hCorr1D_side_1to2p5_Data");
	// double arrContents[] = {1.00, 1.02, 0.95, 1.21, 0.96, 0.93, 0.97, 1.1, 1.1, 1.03};
	double arrContents[] = {0.95, 1.02, 1.02, 1.10, 1.01, 0.96, 0.93, 1.07, 0.98, 0.91};
	double arrErrors[]   = {0.04, 0.03, 0.05, 0.07, 0.03, 0.05, 0.06, 0.12, 0.10, 0.19};
	for (unsigned i = 1; i <= histCorr1D_side_1to2p5_Robin_QCDb->GetNbinsX(); i++) {
		histCorr1D_side_1to2p5_Data->SetBinContent(i,arrContents[i-1]);
		histCorr1D_side_1to2p5_Data->SetBinError(i,arrErrors[i-1]);
	}

	//////////////////////////////////////////
	// Plot individual components together //
	//////////////////////////////////////////
	histCorr1D_side_1to2p5_Robin_QCDb->SetLineWidth(2);
	histCorr1D_side_1to2p5_Robin_QCDb->SetMarkerStyle(21);
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetLineColor(kGreen+2);
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetLineWidth(2);
	// histCorr1D_side_1to2p5_Robin_QCDScatter->SetLineStyle(2);
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetMarkerColor(kGreen+2);
	histCorr1D_side_1to2p5_Robin_QCDScatter->SetMarkerStyle(22);
	histCorr1D_side_1to2p5_Data->SetLineColor(kGreen+2);
	histCorr1D_side_1to2p5_Data->SetLineWidth(2);
	histCorr1D_side_1to2p5_Data->SetMarkerColor(kGreen+2);
	THStack stack("stack","");
	stack.Add(histCorr1D_side_1to2p5_Robin_QCDb);
	// stack.Add(histCorr1D_side_1to2p5_Robin_QCDScatter);
	stack.Add(histCorr1D_side_1to2p5_Data);
	stack.Draw("EPNOSTACK");
	(stack.GetHistogram())->SetXTitle("Bin");
	(stack.GetHistogram())->SetYTitle("Correlation coefficient");
	stack.GetHistogram()->GetXaxis()->SetTitleSize(0.05);
    stack.GetHistogram()->GetXaxis()->SetLabelSize(0.07);
	// stack.GetHistogram()->GetXaxis()->SetTitleOffset(0.05);
    
    stack.GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    stack.GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    // stack.GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
	stack.SetMaximum(2.0);
	stack.SetMinimum(0.0);

	// Add a legend
	// TLegend leg(0.12,0.6,0.4,0.88);
	TLegend leg(0.6,0.67,0.88,0.88);
	leg.SetFillStyle(0);
	leg.AddEntry(histCorr1D_side_1to2p5_Robin_QCDb,"QCD b#bar{b} MC","lp");
	// leg.AddEntry(histCorr1D_side_1to2p5_Robin_QCDScatter,"#splitline{QCD MC (q-g scatter),}{q = b, #bar{b}, c, #bar{c}}","lp");
	leg.AddEntry(histCorr1D_side_1to2p5_Data,"Data","lp");
	leg.Draw();

	// Draw a horizontal line at 1
	double min = histCorr1D_side_1to2p5_Data->GetBinLowEdge(1);
	double max = histCorr1D_side_1to2p5_Data->GetBinLowEdge(histCorr1D_side_1to2p5_Data->GetNbinsX()+1);
	TLine *line = new TLine(min,1,max,1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw();

	TPaveText t(0.15, 0.7, 0.4, 0.8, "NDC");
	t.AddText("Control region A");
    t.SetFillColor(kWhite);
    t.SetBorderSize(0);
    t.Draw();

	// Default canvas has name c1
	c1->SetTicks(1,1); // Put tick marks on top x and right y axes
	c1->SaveAs("../Combined/histCorr1D_side_1to2p5_qcdb_data.pdf");

	THStack* stack_sig = (THStack*) stack.Clone("stack_sig");
	histCorr1D_side_1to2p5_Robin_Signal->SetLineColor(kOrange);
	histCorr1D_side_1to2p5_Robin_Signal->SetLineWidth(2);
	histCorr1D_side_1to2p5_Robin_Signal->SetMarkerColor(kOrange);
	stack_sig->Add(histCorr1D_side_1to2p5_Robin_Signal);
	stack_sig->Draw("EPNOSTACK");
	TLegend* leg_sig = (TLegend*) leg.Clone();
	leg_sig->AddEntry(histCorr1D_side_1to2p5_Robin_Signal, "#splitline{Signal MC}{m_{#varphi} = 8 GeV}","lp");
	leg_sig->Draw();
	line->Draw();
	t.Draw();
	c1->SaveAs("../Combined/histCorr1D_side_1to2p5_signal_qcdb_data.pdf");

	stack.Add(histCorr1D_side_1to2p5_Robin_QCDScatter);
	stack.Draw("EPNOSTACK");
	leg.AddEntry(histCorr1D_side_1to2p5_Robin_QCDScatter,"#splitline{QCD q-g scatter MC}{q = b, #bar{b}, c, #bar{c}}","lp");
	leg.Draw();
	line->Draw();
	t.Draw();
	c1->SaveAs("../Combined/histCorr1D_side_1to2p5.pdf");

	//////////////////////////////////////////////////////////////////
	// Now combine QCD plots - need to reweight for cross-sections //
	//////////////////////////////////////////////////////////////////

	// Need to redo whole correlation coefficient calculation as we have to scale according to # events
	// So get raw #, scale to lumi&xsec, then add, normalise and do as normal.
	TH1D* histM_side_1to2p5_unscaled_QCDb = (TH1D*) fQCDb->Get("hM_side_1to2p5_unnormalised");
	TH1D* histM_side_1to2p5_unscaled_QCDScatter = (TH1D*) fQCDScatter->Get("hM_side_1to2p5_unnormalised");
	std::vector<TH1D*> plots1D;
	plots1D.push_back(histM_side_1to2p5_unscaled_QCDb);
	plots1D.push_back(histM_side_1to2p5_unscaled_QCDScatter);

	TH2D* histM1vsM2_side_1to2p5_unscaled_QCDb = (TH2D*) fQCDb->Get("hM1vsM2_side_1to2p5_unnormalised");
	TH2D* histM1vsM2_side_1to2p5_unscaled_QCDScatter = (TH2D*) fQCDScatter->Get("hM1vsM2_side_1to2p5_unnormalised");
	std::vector<TH2D*> plots2D;
	plots2D.push_back(histM1vsM2_side_1to2p5_unscaled_QCDb);
	plots2D.push_back(histM1vsM2_side_1to2p5_unscaled_QCDScatter);

	///////////////////////////////
	// SET SCALING FACTORS HERE //
	///////////////////////////////
	std::vector<double> scalingFactors;
	scalingFactors.push_back(2.9475); // QCDb
	scalingFactors.push_back(2.6577); // QCDscatter

	// Create combination 2D plot (numerator)
	TH2D* histM1vsM2_side_1to2p5 = (TH2D*) combinePlots(plots2D, scalingFactors);
	drawHistAndSave(histM1vsM2_side_1to2p5, "COLZTEXTE", "histM1vsM2_side_1to2p5_unscaled", directory, app);
	normaliseHist(histM1vsM2_side_1to2p5);

	// Create combination 1D sideband plot (denominator)
	TH1D* histM_side_1to2p5 = (TH1D*) combinePlots(plots1D, scalingFactors);
	drawHistAndSave(histM_side_1to2p5, "HISTE", "histM_side_1to2p5_unscaled", directory, app);
	normaliseHist(histM_side_1to2p5);

	// Create 2D from 1D x 1D
	TH2D* histM1timesM1_side_1to2p5 = create2Dfrom1D(massBins, histM_side_1to2p5);

	// Calculate new correlation coeffs & plot
	TH2D* histM1vsM2_correlations_side_1to2p5 = (TH2D*)histM1vsM2_side_1to2p5->Clone("hM1vsM2_correlations_side_1to2p5");
	histM1vsM2_correlations_side_1to2p5->SetTitle(
		"m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband), (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_side_1to2p5->Divide(histM1timesM1_side_1to2p5);

	// Make 1D plots of unique bins from 2D correlation plot
	TH1D* histCorr1D_side_1to2p5_combo = unique1DBinsFrom2D(histM1vsM2_correlations_side_1to2p5, (nBinsX+1)*nBinsX/2.);

	drawHistAndSave(histM1vsM2_side_1to2p5, "colzTEXTE","M1vsM2_side_1to2p5", directory, app);
	drawHistAndSave(histM1timesM1_side_1to2p5, "colzTEXTE","M1timesM1_side_1to2p5", directory, app);
	drawHistAndSave(histM_side_1to2p5, "HISTE", "M_side_1to2p5", directory, app);
	drawHistAndSave(histM1vsM2_correlations_side_1to2p5, "colzTEXTE","M1vsM2_correlations_side_1to2p5", directory, app);
	histCorr1D_side_1to2p5_combo->SetMaximum(2.0);
	histCorr1D_side_1to2p5_combo->SetMinimum(0);
	drawHistAndSave(histCorr1D_side_1to2p5_combo, "e1", "Correlations1D_side_1to2p5", directory, app);

	// Plot alongside Data
	THStack stack2("stack2","");
	histCorr1D_side_1to2p5_combo->SetLineWidth(2);
	histCorr1D_side_1to2p5_combo->SetMarkerStyle(21);
	stack2.Add(histCorr1D_side_1to2p5_combo);
	stack2.Add(histCorr1D_side_1to2p5_Data);
	histCorr1D_side_1to2p5_Data->SetLineWidth(2);
	stack2.Draw("EPNOSTACK");

	// (stack2.GetHistogram())->SetXTitle("Bin");
	// (stack2.GetHistogram())->SetTitleSize(0.05,"X");
	// (stack2.GetHistogram())->SetLabelSize(0.07,"X");
	// (stack2.GetHistogram())->SetLabelSize(0.05,"Y");
	// (stack2.GetHistogram())->SetYTitle("Correlation coefficient");
	// (stack2.GetHistogram())->SetTitleSize(0.05,"Y");
	// cout << (stack2.GetHistogram())->GetMaximum() << endl;
	// (stack2).SetMaximum(2.0);
	// (stack2).SetMinimum(0);
	setupCorrStuff(stack2.GetHistogram());
	stack2.Draw("EPNOSTACK");

	// Add a legend
	TLegend leg2(0.53,0.63,0.87,0.9);
	// leg2.SetFillColor(kWhite);
	leg2.SetFillStyle(0);
	TH1D blank("","",1,0,1);
	leg2.AddEntry(histCorr1D_side_1to2p5_combo,"Gen. level QCD MC","lp");
	leg2.AddEntry(&blank,"(b#bar{b} + q-g scatter,","");
	leg2.AddEntry(&blank,"q = b, #bar{b}, c, #bar{c})","");
	leg2.AddEntry(histCorr1D_side_1to2p5_Data,"Data","lp");
	leg2.Draw();
	
	line->Draw();

	c1->SetTicks(1,1);
	t.Draw();
	c1->SaveAs("../Combined/histCorr1D_side_1to2p5_combo_allQCD.pdf");
}
