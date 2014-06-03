/**
 * This just post-processes the output from massPlots program, to combine with Alexei's reuslts
 * and other stuff
 */
void massPlotPost() {
	gStyle->SetOptStat(""); // DOES NOTHING AS HIST ALREADY HAS OPT STATS!!!!
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(kWhite);

	// TFile* f = TFile::Open("../QCDScatter_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* fQCDb = TFile::Open("../QCDb_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* fQCDc = TFile::Open("../QCDc_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TFile* fQCDScatter = TFile::Open("../QCDScatter_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT_dR1.root","READ");
	TH1D* hCorr1D_side_1to2p5_Robin_QCDb = (TH1D*) fQCDb->Get("hCorr1D_side_1to2p5");
	TH1D* hCorr1D_side_1to2p5_Robin_QCDc = (TH1D*) fQCDc->Get("hCorr1D_side_1to2p5");
	TH1D* hCorr1D_side_1to2p5_Robin_QCDScatter = (TH1D*) fQCDScatter->Get("hCorr1D_side_1to2p5");
	hCorr1D_side_1to2p5_Robin_QCDb->SetStats(kFALSE);
	hCorr1D_side_1to2p5_Robin_QCDc->SetStats(kFALSE);
	hCorr1D_side_1to2p5_Robin_QCDScatter->SetStats(kFALSE);
	
	// Make Alexei's hist
	TH1D* hCorr1D_side_1to2p5_Alexei = (TH1D*) hCorr1D_side_1to2p5_Robin_QCDb->Clone("hCorr1D_side_1to2p5_Alexei");
	double arrContents[] = {1.00, 1.02, 0.95, 1.21, 0.96, 0.93, 0.97, 1.1, 1.1, 1.03};
	double arrErrors[] = {0.05, 0.03, 0.05, 0.08, 0.04, 0.05, 0.07, 0.13, 0.12, 0.22};
	for (int i = 1; i <= hCorr1D_side_1to2p5_Robin_QCDb->GetNbinsX(); i++) {
		hCorr1D_side_1to2p5_Alexei->SetBinContent(i,arrContents[i-1]);
		hCorr1D_side_1to2p5_Alexei->SetBinError(i,arrErrors[i-1]);
	}

	// Plot 1D hists together
	hCorr1D_side_1to2p5_Robin_QCDScatter->SetLineColor(kGreen+2);
	hCorr1D_side_1to2p5_Robin_QCDScatter->SetMarkerColor(kGreen+2);
	hCorr1D_side_1to2p5_Robin_QCDc->SetLineColor(kOrange+2);
	hCorr1D_side_1to2p5_Robin_QCDc->SetMarkerColor(kOrange+2);
	hCorr1D_side_1to2p5_Alexei->SetLineColor(kRed);
	hCorr1D_side_1to2p5_Alexei->SetMarkerColor(kRed);
	THStack stack("stack","");
	stack.Add(hCorr1D_side_1to2p5_Robin_QCDb);
	// stack.Add(hCorr1D_side_1to2p5_Robin_QCDc);
	stack.Add(hCorr1D_side_1to2p5_Robin_QCDScatter);
	stack.Add(hCorr1D_side_1to2p5_Alexei);
	stack.Draw("EPNOSTACK");
	(stack.GetHistogram())->SetXTitle("Bin");
	(stack.GetHistogram())->SetYTitle("Correlation coefficient");

	TLegend leg(0.65,0.7,0.88,0.88);
	leg.SetFillColor(kWhite);
	leg.AddEntry(hCorr1D_side_1to2p5_Robin_QCDb,"QCDb","lp");
	// leg.AddEntry(hCorr1D_side_1to2p5_Robin_QCDc,"QCDc","lp");
	leg.AddEntry(hCorr1D_side_1to2p5_Robin_QCDScatter,"QCD q-g Scatter","lp");
	leg.AddEntry(hCorr1D_side_1to2p5_Alexei,"Data (Alexei)","lp");
	leg.Draw();

	double min = hCorr1D_side_1to2p5_Alexei->GetBinLowEdge(1);
	double max = hCorr1D_side_1to2p5_Alexei->GetBinLowEdge(hCorr1D_side_1to2p5_Alexei->GetNbinsX()+1);
	TLine *line = new TLine(min,1,max,1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw();

	// Default canvas has name c1
	c1->SetTicks(1,1); // Put tick marks on top x and right y axes
	c1->SaveAs("hCorr1D_side_1to2p5_combined.pdf");
}
