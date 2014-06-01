// // ROOT headers
#include "TStyle.h"

void massPlotPost() {
	gStyle->SetOptStat(""); // DOES NOTHING AS HIST ALREADY HAS OPT STATS!!!!
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(kWhite);

	// TFile* f = TFile::Open("../QCDScatter_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* f = TFile::Open("../QCDb_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TH1D* hCorr1D_side_1to2p5_Robin = (TH1D*) f->Get("hCorr1D_side_1to2p5");
	hCorr1D_side_1to2p5_Robin->SetStats(kFALSE);
	TH1D* hCorr1D_side_1to2p5_Alexei = (TH1D*) hCorr1D_side_1to2p5_Robin->Clone("hCorr1D_side_1to2p5_Alexei");

	double arrContents[] = {1.00, 1.02, 0.95, 1.21, 0.96, 0.93, 0.97, 1.1, 1.1, 1.03};
	double arrErrors[] = {0.05, 0.03, 0.05, 0.08, 0.04, 0.05, 0.07, 0.13, 0.12, 0.22};

	for (int i = 1; i <= hCorr1D_side_1to2p5_Robin->GetNbinsX(); i++) {
		hCorr1D_side_1to2p5_Alexei->SetBinContent(i,arrContents[i-1]);
		hCorr1D_side_1to2p5_Alexei->SetBinError(i,arrErrors[i-1]);
	}

	hCorr1D_side_1to2p5_Robin->Draw("EP");
	hCorr1D_side_1to2p5_Alexei->SetLineColor(kRed);
	hCorr1D_side_1to2p5_Alexei->Draw("EPSAME");

	// Default canvas has name c1
	c1->SaveAs("hCorr1D_side_1to2p5_combined.pdf");
}
